library(tidyverse)
library(rstan)
library(arm)
library(AER) ## for IV
library(lmtest)
library(sandwich)
library(kableExtra)
library(tikzDevice)
library(xtable)
library(texreg)

select <- dplyr::select

source("../regression.r")

### make dataset for analysis
source("code/data.r")

### estimate effects
source("code/estimates.r")

### generate tables and plots
source("code/displayResults.r")
