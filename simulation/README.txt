To replicate* the simulation study in R, follow these steps:

1. Ensure that the following packages are installed in R: (see below for version information)
	- rstan
	- parallel
	- sandwich
	- snow
	- doSNOW
	- utils
	- tidyverse
	- broom
	- kableExtra
	- gridExtra
	- ggpubr
	- tikzDevice

2. Ensure LuaTex is installed

3. Optionally, change the simulation defaults by editing code/runSim.r:
	 - nrep <- 5000  dictates the number of repetitions
	 - ncore <- 10 is for parallel computing--computation is distributed across ncore cores

4. From within this folder (simulation) run the following line of code in R:

> source("code/runSim.r")


* Unfortunately, due to the myriad sources of randomness in the simulation design, we were unable to fix and save all of the randomization seeds. Therefore, replication results will not be exactly the same as what was reported in the paper.


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
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tidyr_1.3.1         doSNOW_1.0.20       iterators_1.0.14    foreach_1.5.2      
 [5] snow_0.4-4          purrr_1.0.4         sandwich_3.1-1      rstan_2.32.6       
 [9] StanHeaders_2.32.10 dplyr_1.1.4        

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       compiler_4.5.0     tidyselect_1.2.1   Rcpp_1.0.13-1      gridExtra_2.3     
 [6] scales_1.4.0       lattice_0.22-6     ggplot2_4.0.0      R6_2.6.1           generics_0.1.3    
[11] tibble_3.2.1       pillar_1.10.1      RColorBrewer_1.1-3 rlang_1.1.5        inline_0.3.20     
[16] S7_0.2.0           RcppParallel_5.1.9 cli_3.6.4          magrittr_2.0.3     grid_4.5.0        
[21] lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0         farver_2.1.2       QuickJSR_1.4.0    
[26] codetools_0.2-20   zoo_1.8-14         stats4_4.5.0       pkgbuild_1.4.5     colorspace_2.1-1  
[31] matrixStats_1.4.1  tools_4.5.0        loo_2.8.0          pkgconfig_2.0.3   
