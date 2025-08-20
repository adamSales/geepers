# Replication material for
# "GEEPERs: Principal Stratification using Principal Scores and Stacked Estimating
# Equations"


Replications will not be exact, mostly because of randomness in the simulation and some of the estimation, and also because of very slight differences (e.g. one additional subject) between the publicly available data and the dataset used in the paper.


## Replication the Simulation
To replicate the simulation in R, install the required packages (see below) then (from the main folder) run: (this will take a loonnnng time)

In R:

```
> source("code/simulation/runSim.r")
```
from the command line:

```
> Rscript code/simulation/runSim.r
```


## Replicating the OPTS (Periodontal Disease) Analysis
To replicate the OPTS analysis, run:

```
> source("code/OPTS/opts.r")
```

## Replicating the Bottom-Outer Study

### Obtaining the Data
To obtain the data, follow the instructions at (https://osf.io/skghn).
After doing so, you will receive a link to a new OSF link. 
From that OSF site, download the following files:

- Files/Data/assessments/assess_student.csv
- Files/Data/ASSISTments/assist_student_problem.csv
- Files/Data/student/student_attendance.csv
- Files/Data/student/student_demo.csv
- Files/Data/student/student_roster.csv

Place them all in the "data" folder. 

Then, from this folder, run:

```
> source("code/fh2t/replicate.r")
```

















```
R version 4.5.0 (2025-04-11)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.5

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggpubr_0.6.0         gridExtra_2.3        broom_1.0.7          lubridate_1.9.4      stringr_1.5.1       
 [6] tidyverse_2.0.0      tikzDevice_0.12.6    coefplot_1.2.8       tibble_3.2.1         purrr_1.0.4         
[11] xtable_1.8-4         tableone_0.13.2      texreg_1.39.4        kableExtra_1.4.0     rstan_2.32.6        
[16] StanHeaders_2.32.10  estimatr_1.0.4       missForest_1.5       randomForest_4.7-1.2 arm_1.14-4          
[21] lme4_1.1-35.5        Matrix_1.7-3         MASS_7.3-65          forcats_1.0.0        lmtest_0.9-40       
[26] zoo_1.8-14           sandwich_3.1-1       tidyr_1.3.1          readr_2.1.5          ggplot2_3.5.1       
[31] dplyr_1.1.4         

loaded via a namespace (and not attached):
 [1] DBI_1.2.3          inline_0.3.20      rlang_1.1.5        magrittr_2.0.3     matrixStats_1.4.1 
 [6] compiler_4.5.0     mgcv_1.9-1         loo_2.8.0          systemfonts_1.1.0  vctrs_0.6.5       
[11] reshape2_1.4.4     crayon_1.5.3       pkgconfig_2.0.3    fastmap_1.2.0      backports_1.5.0   
[16] labeling_0.4.3     utf8_1.2.4         rmarkdown_2.29     tzdb_0.4.0         nloptr_2.1.1      
[21] itertools_0.1-3    ragg_1.3.3         xfun_0.52          parallel_4.5.0     R6_2.6.1          
[26] stringi_1.8.4      car_3.1-3          boot_1.3-31        Rcpp_1.0.13-1      iterators_1.0.14  
[31] knitr_1.49         filehash_2.4-6     timechange_0.3.0   tidyselect_1.2.1   rstudioapi_0.17.1 
[36] abind_1.4-8        codetools_0.2-20   pkgbuild_1.4.5     doRNG_1.8.6        lattice_0.22-6    
[41] plyr_1.8.9         withr_3.0.2        coda_0.19-4.1      evaluate_1.0.1     survival_3.8-3    
[46] RcppParallel_5.1.9 survey_4.4-2       xml2_1.3.6         useful_1.2.6.1     pillar_1.10.1     
[51] carData_3.0-5      rngtools_1.5.2     foreach_1.5.2      stats4_4.5.0       generics_0.1.3    
[56] hms_1.1.3          munsell_0.5.1      scales_1.3.0       minqa_1.2.8        glue_1.8.0        
[61] tools_4.5.0        ggsignif_0.6.4     cowplot_1.1.3      grid_4.5.0         mitools_2.4       
[66] QuickJSR_1.4.0     colorspace_2.1-1   nlme_3.1-168       Formula_1.2-5      cli_3.6.4         
[71] textshaping_0.4.1  viridisLite_0.4.2  svglite_2.1.3      gtable_0.3.6       rstatix_0.7.2     
[76] digest_0.6.37      farver_2.1.2       htmltools_0.5.8.1  lifecycle_1.0.4    httr_1.4.7        
```
