
data{
 int<lower=1> nc; // number of control subjects
 int<lower=1> nt; // number of treated subjects
 int<lower=0> ncovU; // number of covariates
 int<lower=0> ncovY; // number of covariates  

 real YctlY[nc]; // control outcomes
 real YtrtY[nt]; // treatment outcomes
 
 matrix[nc,ncovU] XctlU; // covariate matrix for controls
 matrix[nt,ncovU] XtrtU; //covariate matrix for treated

 matrix[nc,ncovY] XctlY; // covariate matrix for controls
 matrix[nt,ncovY] XtrtY; //covariate matrix for treated

 int<lower=0,upper=1> complier[nt]; // 
}


parameters{
 real alphaTBO; // Y-intercept for treatment Bottom Outers
 real alphaTNBO; // Y-intercept for treatment Not  Bottom Outers
 real alphaCBO; // Y-intercept for control  Bottom Outers
 real alphaCNBO; // Y-intercept for control Not Bottom Outers

 real alphaU; // intercept for usage model
 vector[ncovU] betaU; // coefficients for covariates in usage model
 vector[ncovY] betaY; // coefficients for covariates in outcome model 

 // residual standard deviation:
 real<lower=0> sigTBO;
 real<lower=0> sigTNBO;
 real<lower=0> sigCBO;
 real<lower=0> sigCNBO;
}


transformed parameters{
 //declare new parameters:
 real complierATE; //Avg. Effect for skippers
 real neverTakerATE; //Avg. Effect for firsters
 real ATEdiff; //Difference btw Avg Effects
 vector[nc] piC; //Pr(complier) for controls
 vector[nt] piT; //Pr(complier) for treateds

 //define new parameters:
 complierATE=alphaTBO-alphaCBO; // effect of treatement on bottom outers
 neverTakerATE=alphaTNBO-alphaCNBO; // effect of treatment on not bottom outers
 ATEdiff=complierATE-neverTakerATE; // differences in effects
 piC=inv_logit(alphaU+XctlU*betaU);
 piT=inv_logit(alphaU+XtrtU*betaU);// prop of bing a bottom outer
}


model{
 // vectors of intercepts and residual SDs for treated
 // students. useful for vectorizing:
 vector[nt] alphaT;
 vector[nt] sigT;

 // The ? functions as if-else A?B:C returns B if A=1 and C otherwise https://mc-stan.org/docs/2_29/reference-manual/conditional-statements.html
 for(i in 1:nt){
   alphaT[i] = complier[i]?alphaTBO:alphaTNBO;
   sigT[i] = complier[i]?sigTBO:sigTNBO;
 }

// std_normal() prior for everything?
alphaTBO~std_normal();
alphaTNBO~std_normal();
alphaCBO~std_normal();
alphaCNBO~std_normal();
alphaU~std_normal();
betaU~std_normal();
betaY~std_normal();
sigTBO~std_normal();
sigTNBO~std_normal();
sigCBO~std_normal();
sigCNBO~std_normal();



 complier~bernoulli(piT); // model for who is a bottom outer
 YtrtY~normal(alphaT+XtrtY*betaY,sigT); // model for outcomes in treatment group


for(i in 1:nc)
 target += log_sum_exp(
  log(piC[i]) + normal_lpdf(YctlY[i]| alphaCBO+XctlY[i,]*betaY,sigTBO),
  log((1-piC[i])) + normal_lpdf(YctlY[i] |alphaCNBO+XctlY[i,]*betaY,sigTNBO));
} 