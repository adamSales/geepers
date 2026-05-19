To replicate the video scaffolding example in R, follow these steps:

1. From https://osf.io/59shv/files/osfstorage, download:
	 - experiment_dataset_2021-09-23.zip (do not unzip this)
	 - experiment_descriptions.csv
	 - processed_experiment_data.csv
	 and place them in the "rawdata" folder

2. Ensure that the following packages are installed in R: (see below for version information)
	 - tidyverse
	 - rstan
	 - arm
	 - AER 
	 - lmtest
	 - sandwich
	 - kableExtra
	 - tikzDevice
	 - xtable
	 - texreg

3. Ensure LuaTex is installed

4. From within this folder (videoScaffolding) run the following line of code in R:

> source("code/replicate.r")



R session information

> sessionInfo()
R version 4.5.0 (2025-04-11)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.7.4

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] texreg_1.39.4       xtable_1.8-4        tikzDevice_0.12.6   kableExtra_1.4.0    AER_1.2-14         
 [6] survival_3.8-3      sandwich_3.1-1      lmtest_0.9-40       zoo_1.8-14          car_3.1-3          
[11] carData_3.0-5       arm_1.14-4          lme4_1.1-35.5       Matrix_1.7-3        MASS_7.3-65        
[16] rstan_2.32.6        StanHeaders_2.32.10 lubridate_1.9.4     forcats_1.0.0       stringr_1.5.1      
[21] dplyr_1.1.4         purrr_1.0.4         readr_2.1.5         tidyr_1.3.1         tibble_3.2.1       
[26] ggplot2_4.0.0       tidyverse_2.0.0    

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1   filehash_2.4-6     viridisLite_0.4.2  farver_2.1.2       loo_2.8.0         
 [6] S7_0.2.0           fastmap_1.2.0      digest_0.6.37      timechange_0.3.0   lifecycle_1.0.4   
[11] processx_3.8.4     magrittr_2.0.3     compiler_4.5.0     rlang_1.1.5        tools_4.5.0       
[16] knitr_1.49         labeling_0.4.3     bit_4.5.0.1        pkgbuild_1.4.5     xml2_1.3.6        
[21] RColorBrewer_1.1-3 abind_1.4-8        withr_3.0.2        grid_4.5.0         stats4_4.5.0      
[26] colorspace_2.1-1   inline_0.3.20      scales_1.4.0       cli_3.6.4          crayon_1.5.3      
[31] rmarkdown_2.29     ragg_1.3.3         generics_0.1.3     RcppParallel_5.1.9 rstudioapi_0.17.1 
[36] httr_1.4.7         tzdb_0.4.0         minqa_1.2.8        splines_4.5.0      parallel_4.5.0    
[41] matrixStats_1.4.1  vctrs_0.6.5        boot_1.3-31        callr_3.7.6        hms_1.1.3         
[46] bit64_4.5.2        Formula_1.2-5      systemfonts_1.1.0  glue_1.8.0         nloptr_2.1.1      
[51] ps_1.8.1           codetools_0.2-20   stringi_1.8.4      gtable_0.3.6       QuickJSR_1.4.0    
[56] pillar_1.10.1      htmltools_0.5.8.1  R6_2.6.1           textshaping_0.4.1  vroom_1.6.5       
[61] evaluate_1.0.1     lattice_0.22-6     Rcpp_1.0.13-1      svglite_2.1.3      coda_0.19-4.1     
[66] gridExtra_2.3      nlme_3.1-168       xfun_0.52          pkgconfig_2.0.3   
