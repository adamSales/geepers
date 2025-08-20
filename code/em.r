## adapted from chatGPT

## data includes Y, Z, x1,x2, and S[Z==1]

EM_gaussian_mixture <- function(y, p, mu01_init, mu00_init, sigma2_init = 1, tol = 1e-6, max_iter = 1000) {
  n <- length(y)
  mu01 <- mu01_init
  mu00 <- mu00_init
  sigma2 <- sigma2_init

  loglik <- function(mu01, mu00, sigma2) {
    comp1 <- dnorm(y, mean = mu01, sd = sqrt(sigma2))
    comp2 <- dnorm(y, mean = mu00, sd = sqrt(sigma2))
    sum(log(p * comp1 + (1 - p) * comp2))
  }

  ll_tot <- numeric(max_iter+1)
  ll_tot[1] <- ll_old <- loglik(mu01, mu00, sigma2)

  for (iter in 1:max_iter) {
    # E-step
    comp1 <- dnorm(y, mean = mu01, sd = sqrt(sigma2))
    comp2 <- dnorm(y, mean = mu00, sd = sqrt(sigma2))
    tau <- (p * comp1) / (p * comp1 + (1 - p) * comp2)

    # M-step
    mu01 <- sum(tau * y) / sum(tau)
    mu00 <- sum((1 - tau) * y) / sum(1 - tau)
    sigma2 <- sum(tau * (y - mu01)^2 + (1 - tau) * (y - mu00)^2) / n

    # Check convergence
    ll_tot[iter+1] <- ll_new <- loglik(mu01, mu00, sigma2)
    if (abs(ll_new - ll_old) < tol) {
      break
    }
    ll_old <- ll_new
  }


  list(mu01 = mu01, mu00 = mu00, sigma2 = sigma2,
       tau = tau, loglik = ll_new, ll_tot=ll_tot[1:(iter+1)],iter = iter)
}

PS_EM <- function(data,psMod=NULL,tol = 1e-6, max_iter = 1000) {

    ## estimate principal scores
    if(is.null(psMod)) psMod=glm(S~x1+x2,data=data,subset=Z==1,family=binomial)

    ## outcome model for treatment group
    out1 <- lm(Y~S+x1+x2,data=data, subset=Z==1)

    mu11 <- coef(out1)["(Intercept)"]+coef(out1)["S"]
    mu10 <- coef(out1)["(Intercept)"]

    ## estimate intercepts for control group
    dat0 <- data[data$Z==0,]
    ctl_sol <- EM_gaussian_mixture(
        y= with(dat0, Y-coef(out1)["x1"]*x1-coef(out1)["x2"]*x2),
        p= predict(psMod,newdata=dat0,type="response"),
        mu01_init=mu11,
        mu00_init=mu10,
        sigma2_init=summary(out1)$sigma^2,
        tol=tol,
        max_iter=max_iter)

    out <- c(eff0=unname(mu10-ctl_sol$mu00),
             eff1=unname(mu11-ctl_sol$mu01))
    attr(out,"psMod") <- psMod
    attr(out,"out1") <- out1
    attr(out,"ctl_sol") <- ctl_sol
    out
}

