#####################################################################################################
# Fuzzy Regression Discontinuity Designs with Multiple Control Groups Under One-Sided Noncompliance: 
# Evaluating Extended Time Accommodations
#
# Youmi Suk & Yongnam Kim
#####################################################################################################

# load packages
library(MatchIt)
library(rdrobust)

# data generating models
DGP <- function(n=2000, beta=2, error.sd=0.01, CausalAssumptionsViolation=NULL) {
  
  repeat {
    W <- runif(n, -1, 1) # assuming W is measured.
    X <- runif(n, min=-0.4, max=1.4) # running variable
    U <- runif(n, -1, 1) # assuming U is unmeasured.
    
    A <- as.numeric(I(X <=0)) # eligibility
    
    Z.prob <- 1/20*(5 + 6*(A) - 0.5*X  + 1*W - 3*U) + rnorm(n, 0, error.sd)
    Z <- as.numeric(Z.prob > 0.5)
    
    Trt.prob <- 1/20*(7 + 3*(A) - 4*(Z==0) - 0.5*X  + 1*W  + 1*U) + rnorm(n, 0, error.sd)
    Trt <- as.numeric(Trt.prob > 0.5)
    
    if (sum(table(A, Trt)[[1, 2]] == 0 & table(Z, Trt)[[1, 2]] == 0) == 1) {
      break
    }
    
  }   
  
  if (is.null(CausalAssumptionsViolation)) {
    Y0 <- 8 + 3*X + 1*W + 1*U + 0.8*as.numeric(I(X^2)) + rnorm(n)
    Y1 <- Y0 + 1 - 0.5*X
  } else if (CausalAssumptionsViolation == "Exclusion") {
    Y0 <- 8 + 3*X + 1*W + 1*U + 0.8*as.numeric(I(X^2)) + rnorm(n) + 0.5*A   
    Y1 <- Y0 + 1 - 0.5*X
  }
  
  Y <- ifelse(Trt==1, Y1, Y0)
  
  dat <- data.frame(A, Z, W, Trt, X, Trt.prob, Y1, Y0, Y, U)
  
  return(dat)
}

set.seed(1)
reps = 500 # number of replications
rlst_list <- rlst_bw <- list()

for (i in 1:2) {
  rlst_list[[i]] <- matrix(NA, nrow=reps, ncol=9)
}

for (c in 1:2) {
  
  if (c == 1) {
    violation = NULL
  } else if (c == 2) {
    violation = "Exclusion"
  } 
  
  for (i in 1:reps) {
    
    print(paste("Condition: ", c, " & Replication ", i))  
    
    dat <- DGP(n = 10000, CausalAssumptionsViolation=violation) 
    
    # --- Fuzzy RD Design
    bw <- rdbwselect(y = dat$Y, x = dat$X, fuzzy = dat$Trt, kernel = "tri", p = 2)
    rd_bw <- bw$bws[1]   
    
    t1 <- lm(Trt ~ A + X + I(X^2)  + A:X,  dat) # IV - correct specification
    y1 <- lm(Y ~ A + X + I(X^2)  + A:X, dat)
    
    t1_2 <- lm(Trt ~ A + X  + A:X,  dat) # IV - incorrect specification
    y1_2 <- lm(Y ~ A + X + A:X, dat)
    
    # --- subset data
    dat1 <- dat[dat$Z==1, ]
    dat2 <- dat[dat$Z==dat$Trt, ] 
    dat1.X <- dat[dat$Z==1 & dat$X >= -rd_bw & dat$X <= 0,] 
    dat2.X <- dat[dat$Z==dat$Trt & dat$X >= -rd_bw  & dat$X <= 0, ]
    dat1.trtX <- dat1[dat1$Trt==1 & dat1$X >= -rd_bw  & dat1$X <= 0, ]
    dat2.trtX <- dat2[dat2$Trt==1 & dat2$X >= -rd_bw  & dat2$X <= 0, ]
    
    # --- Outcome regression
    out1_cnt <- lm(Y ~ W, data=dat1, subset= Trt==0 & X <= 0 & X >= -rd_bw)
    out2_cnt <- lm(Y ~ W, data=dat2, subset= Trt==0 & X <= 0 & X >= -rd_bw)
    
    out1.Y0hat <- predict(out1_cnt, data.frame(W=dat1.trtX$W, X=dat1.trtX$X)) 
    out2.Y0hat <- predict(out2_cnt, data.frame(W=dat2.trtX$W, X=dat1.trtX$X))
    
    out1.ATT <- mean(dat1.trtX$Y - out1.Y0hat) 
    out2.ATT <- mean(dat2.trtX$Y - out2.Y0hat)
    
    # --- PS Matching
    mF1.X <- matchit(Trt ~ W, data = dat1.X, replace=TRUE,
                     estimand = "ATT", distance = "mahalanobis")  
    mF2.X <- matchit(Trt ~ W, data = dat2.X, replace=TRUE,
                     estimand = "ATT", distance = "mahalanobis") 
    
    md1.X <- match.data(mF1.X)
    md2.X <- match.data(mF2.X)
    
    mat1X.ATT <- coef(lm(Y ~ Trt, data = md1.X, weights = weights))[2]
    mat2X.ATT <- coef(lm(Y ~ Trt, data = md2.X, weights = weights))[2]
    
    # --- PS Weighting
    ps1_X <- glm(Trt ~ W, data = dat1.X, family = 'binomial')
    ps2_X <- glm(Trt ~ W, data = dat2.X, family = 'binomial')
    
    dat1.X$ps <- predict(ps1_X, type="response")  
    dat2.X$ps <- predict(ps2_X, type="response")  
    
    dat1.X$ipwt <- with(dat1.X, Trt  + (1 - Trt) * ps / (1 - ps)) # IPW weight for the upper bound
    dat2.X$ipwt <- with(dat2.X, Trt  + (1 - Trt) * ps / (1 - ps)) # IPW weight for the lower bound
    
    ipwt1X.ATT <- coef(lm(Y ~ Trt, data = dat1.X, weights = ipwt))[2] 
    ipwt2X.ATT <- coef(lm(Y ~ Trt, data = dat2.X, weights = ipwt))[2]
    
    temp_rlst <- as.numeric(c(nrow(dat[dat$Trt == 1 & dat$X >= -rd_bw & dat$X <= 0, ]), # sample size near the cutoff
                              c(y1$coefficients[2]/t1$coefficients[2], y1_2$coefficients[2]/t1_2$coefficients[2]), # fuzzy RD estimate from correct specification, and one from incorrect specification

                              out1.ATT, # upper bound based on outcome regression
                              mat1X.ATT, ipwt1X.ATT, # upper bound based on matching and upper bound based on weighting                    
                              
                              out2.ATT, # lower bound based on outcome regression 
                              mat2X.ATT, ipwt2X.ATT)) # lower bound based on matching and lower bound based on weighting 
    
    
    rlst_list[[c]][i,] <- temp_rlst
    
  }
}  

rlst1_dat <- rlst_list[[1]] # no violation of causal assumptions
rlst2_dat <- rlst_list[[2]] # violation of the exclusion restriction assumption