# Fuzzy Regression Discontinuity Designs with Multiple Control Groups Under One-Sided Noncompliance: Evaluating Extended Time Accommodations

Youmi Suk<sup>1</sup> and Yongnam Kim<sup>2</sup>

<sup>1</sup> Teachers College Columbia University  
<sup>2</sup> Seoul National University


## Overview

Since its conception by Thistlethwaite and Campbell (1960), a regression discontinuity (RD) design has continuously evolved and progressed in various aspects. This paper explores a novel theoretical advancement in fuzzy RD designs by combining them with multiple control groups under one-sided noncompliance. Within this approach, we estimate the average treatment effect on the treated at the cutoff from the fuzzy RD design and construct its bounds derived from multiple control groups. Using the bounds as a sensitivity check, we examine whether the underlying causal or statistical assumptions for the fuzzy RD design are warranted. We verify the effectiveness of our approach via a simulation study and demonstrate its application by studying the effect of extended time accommodations using data from the National Assessment of Educational Progress. 

For more details, see [our paper](https://doi.org/10.31234/osf.io/sa96g). Here, we provide `R` codes to reproduce our simulation study. 

## Simulation Study

* `SimuCode.R`  

   This `R` file includes a function for generating simulated data (i.e., `DGP`) where `CausalAssumptionsViolation` controls the violation of a causal assumption, specifically the exclusion restriction assumption. Additionaly, this file contains simulation codes for different conditions. For more information, see [our paper](https://doi.org/10.31234/osf.io/sa96g).
