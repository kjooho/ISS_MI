# ISS_MI

R code for "Scalable and Robust Multiple Imputation for Case-Cohort Studies via Influence Function-Based Supersampling" by Jooho Kim, Takumi Saegusa, and Yei Eun Shin.

## Directory Structure

### **R**
This directory contains R functions for implementing multiple imputation for supersampled case-cohort and stratified supersampled case-cohort designs using influence function-based supersampling (ISS).

- **supersampling.R**  
  R function for random supersampling (RSS) and influence function-based supersampling (ISS).

- **addph2var.R**  
  R function for phase-2 variance estimation in case-cohort studies, implementing Lin–Ying, Borgan II, and Samuelsen (2007) variance estimators.

- **smcfcs.cc.R**  
  A modified version of the substantive model compatible FCS algorithm from the `smcfcs` package (Bartlett et al., 2015), adapted for case-cohort and supersampled case-cohort designs.

- **generate.data.R**  
  R function for simulating full cohort datasets for demonstration.

### **demo/**
- **casecohort_demo.R**  
  A demonstration script showing how to perform multiple imputation and Cox proportional hazards analysis in supersampled case-cohort studies using ISS.
