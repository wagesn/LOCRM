# Local Continumal Reassessment Methods

R codes to implement local continual reassessment methods for dose finding and optimization in drug-combination trials.

# Description

Due to the limited sample size and large dose exploration space, obtaining a desirable dose combination is a challenging task in the early development of combination treatments for cancer patients. Most existing designs for optimizing the dose combination are model-based, requiring significant efforts to elicit parameters or prior distributions. Model-based designs also rely on intensive model calibration and may yield unstable performance in the case of model misspecification or sparse data. We propose to employ local, underparameterized models for dose exploration to reduce the hurdle of model calibration and enhance the design robustness. Building upon the framework of the partial ordering continual reassessment method (POCRM), we develop local data-based CRM (LOCRM) designs for identifying the maximum tolerated dose combination (MTDC), using toxicity only, and the optimal biological dose combination (OBDC), using both toxicity and efficacy, respectively. The LOCRM designs only model the local data from neighboring dose combinations. Therefore, they are flexible in estimating the local space and circumventing unstable characterization of the entire dose-exploration surface. Our simulation studies show that our approach has competitive performance compared to widely used methods for finding MTDC, and it has advantages over existing model-based methods for optimizing OBDC.

# Functions

The repository includes two functions:

- locrm.R: The R code that includes the function ```locrm()``` to obtain the operating characteristics of the LOCRM design for single agent trials with late-onset toxicities by simulating trials.
  
  ```rscript
  locrm(p.true.tox,target.tox,cutoff.tox,ntrial,nmax,cohortsize)
  ```
  
- locrm12.R: The R code that includes the function ```locrm12()``` to obtain the operating characteristics of the time-to-event model-assisted design for single agent trials with late-onset toxicities by simulating trials.
  
  ```rscipt
  locrm12(p.true.tox,p.true.eff,target.tox,cutoff.tox,target.eff,cutoff.eff,ntrial,nmax,cohortsize)
  ```
  

# Inputs

- `p.true.tox`: A matrix containing the true toxicity probabilities of the investigational dose combinations.
  
- `p.true.eff`: A matrix containing the true efficacy probabilities of the investigational dose combinations.
  
- `target.tox`: The upper limit of toxicity rate, e.g., `target.tox<-0.3`.
  
- `target.eff`: The lower limit of efficacy rate, e.g., `target.eff<-0.2`.
  
- `cutoff.tox`: The cutoff to eliminate an overly toxic dose for safety.
  
- `cutoff.eff`: The cutoff to eliminate an ineffective dose for futility.
  
- `ntrial`: The total number of trial to be simulated.
  
- `nmax`: The maximum sample size.
  
- `cohortsize`: The cohort size.
  

# Outputs

- `locrm()` will return the operating characteristics of the LOCRM design as a data frame, including:

  ```
   (1) selection percentage at each dose combination (sel);  
   (2) the number of toxicities treated at each dose combination (tox);  
   (3) the number of patients treated at each dose combination (pts);  
   (4) the average number of toxicities (ntox);  
   (5) the average number of patients (npts);  
   (6) the percentage of early stopping without selecting the MTDC (p.stop)  
   (7) the average DLT rate (dlt.rate)
  ```
  
- `locrm12()` will return the operating characteristics of the LOCRM12 design as a data frame, including:
  
  ```
   (1) selection percentage at each dose combination (sel);  
   (2) the number of toxicities treated at each dose combination (tox);  
   (3) the number of efficacies treated at each dose combination (eff);  
   (4) the number of patients treated at each dose combination (pts);  
   (5) the average number of toxicities (ntox);  
   (6) the average number of efficacies (neff);  
   (7) the average number of patients (npts);  
   (8) the percentage of early stopping without selecting the OBDC (p.stop)
  ```

# Example

We consider use the LOCRM design as an illustration.

- Suppose the upper limit of toxicity rate is 35%, the lower limit of efficacy rate is 20%. A maximum of 51 patients will be recruited in cohorts of 3. Suppose 5 dose levels of agent A and 3 dose levels of agent B are considered are considered. We can use the following code to simulate the scenario 1 in Table 3.
  
```rscript
  p.true.tox <- matrix(c(0.05,0.15,0.30,0.45,0.55,  0.15,0.30,0.45,0.55,0.65,  0.30,0.45,0.55,0.65,0.75),nrow = 3, byrow = TRUE)
  p.true.eff <- matrix(c(0.05,0.25,0.50,0.55,0.60,  0.25,0.50,0.55,0.60,0.65,  0.50,0.55,0.60,0.65,0.70),nrow = 3, byrow = TRUE)
  locrm12(p.true.tox,p.true.eff,target.tox=0.35,cutoff.tox=0.85,target.eff=0.2,cutoff.eff=0.9,nmax=51,cohortsize=3,ntrial=5000)
  
  -----------------------output------------------------

$sel
      [,1]  [,2]  [,3] [,4] [,5]
[1,]  0.92  3.84 18.38 4.68 0.26
[2,]  6.50 24.22  9.42 0.44 0.00
[3,] 20.56  8.88  0.88 0.02 0.00

$tox
       [,1]   [,2]   [,3]   [,4]   [,5]
[1,] 0.1938 0.5586 2.0422 1.0324 0.1054
[2,] 0.7246 2.8496 1.9618 0.3190 0.0256
[3,] 2.3824 2.4500 0.6954 0.1212 0.0072

$eff
       [,1]   [,2]   [,3]   [,4]   [,5]
[1,] 0.1938 0.5586 2.0422 1.0324 0.1054
[2,] 0.7246 2.8496 1.9618 0.3190 0.0256
[3,] 2.3824 2.4500 0.6954 0.1212 0.0072

$pts
       [,1]   [,2]   [,3]   [,4]   [,5]
[1,] 3.6834 3.6726 6.7446 2.2950 0.1968
[2,] 4.8282 9.4038 4.3812 0.5850 0.0402
[3,] 7.8864 5.3640 1.2738 0.1896 0.0096

$ntox
[1] 15.469

$neff
[1] 22.364

$npts
[1] 50.554

$p.stop
[1] 1
```

  # Authors and Reference
  * Jingyi Zhang, Fangrong Yan, Nolan A. Wages, and Ruitao Lin
  * Zhang, J., Yan, F., Wages, N. A., and Lin, R. (2023) “Local Continual Reassessment Methods for Dose Finding and Optimization in Drug-Combination Trials”, Statistical Methods in Medical Research, revision submitted to journal.

