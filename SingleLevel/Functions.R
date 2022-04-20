#Functions.R for single level analyses

##################### Packages - Starts #############################
library(ggmcmc)
#library(R2jags)
library(robustHD)
library(ggpubr)
library(MASS)
#library(psych)
#library(xlsx)
library(rstan)
library(loo)
options(mc.cores = 4)
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
##################### Packages - Ends #############################

#Functions
my.bin <- function(obs.data, syn.data){
  #This function throws continuous PARED into bins. Bins are the same with original pisa cycle
  values <- as.numeric(levels(factor(obs.data))) #check what are the values in the original pisa cycle
  break.point <- rep(0, length(values)) #generate a break point vector
  for(i in 1:(length(values))){ #prepare the break points 
    break.point[i] <- (values[i]+values[i+1])/2} #each break point is the mean of it self and the next value
  break.point[is.na(break.point)] <- Inf #Inf, because we don't know the maximum value produced by mvnorm function. 
  break.point <- c(0, break.point) #add initial value as 0
  data <- cut(syn.data, breaks=break.point, include.lowest=TRUE, right=FALSE, labels=values) #cut is a base function
  return(as.numeric(as.character(data))) #save from being a factor and turn into numeric values
}

generate.pisa <- function(pisa.data, n=n, cont.var = c(2,3), res.std.err = 84, coef) {
  ## This function generates PISA cycles with "Female" "LangAtHome" "PARED" "CULTPOSS" "HEDRES" "HOMEPOS" "PV1MATH" by taking into account of their correlations and their structure
  Sig_s <-  cov(pisa.data[,cont.var]) #cont.var = continuous variables - we produce PARED as continuous and we'll turn into bins
  S1<-mvrnorm(n,c(colMeans(pisa.data)[cont.var]),Sig_s) #mean values are coming from the original dataset
  Female<-rbinom(n,1,mean(pisa.data$Female))#binary variable probabilities coming from the original dataset percentages
  IMMIG<-rbinom(n,1,mean(pisa.data$IMMIG))#binary variable probabilities coming from the original dataset percentages
  S1[,"PARED"] <- my.bin(obs.data=pisa.data$PARED, syn.data=S1[,"PARED"]) #to bin PARED, we are using our own function
  resi<-rnorm(n,0,res.std.err) #add some residual standard error to PV1MATH
  syn.dat <- as.matrix(cbind(1, Female, S1, IMMIG)) #save the data as matrix and 1 for intercept - unfactor function comes from varhandle package
  PV1MATH <- syn.dat %*% coef + resi #Generate PV!MATH based on our coefficients, I found 84 from lm function
  final.syn.dat <- data.frame(Female, S1, IMMIG, PV1MATH) #generated dataset with PV1MATH
  return(final.syn.dat)
}

#measurement

MSE <- function(y.true, y.pred){(y.true-y.pred)^2}
bias <- function(y.true, y.pred){(y.pred-y.true)}
abs.bias <- function(y.true, y.pred){abs(y.pred-y.true)}
per.bias <- function(y.true, y.pred){(y.pred-y.true)/y.true}
abs.per.bias <- function(y.true, y.pred){abs((y.pred-y.true)/y.true)}

data.stan <- function(dat_combined) {
  dat <- dat_combined[rank(dat_combined$cycle, ties.method = 'first'), ]
  x <- cbind(Intercept = 1, dat[, 1:4])
  y <- dat$PV1MATH
  cycle <- dat$cycle
  N <- nrow(x)
  M <- ncol(x)
  C <- length(unique(cycle))
  NN <- sum(cycle == 1)
  list(N = N, M = M, C = C, NN = NN, x = x, y = y, cycle = cycle)
}

##### non inf

modelstring <- '
data {
  int<lower = 1> N;                    // number of students
  int<lower = 1> M;                    // number of covariates
  int<lower = 1> C;                    // number of cycles
  int<lower = 1> NN;                   // number of students in the current cycle
  matrix[N, M] x;                      // covariates
  vector[N] y;                         // outcomes
  int<lower = 1, upper = C> cycle[N];  // cycles
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;
  
  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;
  
  // x[, 1] is the intercept
  x_std[, 1] = x[, 1];
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
}

parameters {
  vector[M] beta_std;
  real<lower = 0> sigma_y_std;
}

model {
  beta_std ~ normal(0, 100);
  sigma_y_std ~ cauchy(0, 2);
  
  y_std ~ normal(x_std * beta_std, sigma_y_std);
}

generated quantities {
  vector[M] beta;
  real<lower = 0> sigma_y = sigma_y_std * sd_y;
  vector[NN] log_lik;
  
  beta[1] = sd_y * beta_std[1] + mu_y;
  for (m in 2:M) {
    beta[m] = sd_y / sd_x[m] * beta_std[m];
    beta[1] -= beta[m] * mu_x[m];
  }
  
  for (n in 1:NN)
    log_lik[n] = normal_lpdf(y[n] | x[n, ] * beta, sigma_y);
}
'

my.blr <- function(dat_combined, n.chains = 4) {
  writeLines(modelstring, con = 'modelBLR.stan')
  data <- data.stan(subset(dat_combined, cycle == 1))
  stan('modelBLR.stan', data = data, iter = 20000, chains = 4)
}

my.blr.pool <- function(dat_combined, n.chains = 4) {
  writeLines(modelstring, con = 'modelBLR.stan')
  data <- data.stan(dat_combined)
  stan('modelBLR.stan', data = data, iter = 20000, chains = 4)
}

######## Informative Model
modelstring.inf <- '
data {
  int<lower = 1> N;                    // number of students
  int<lower = 1> M;                    // number of covariates
  int<lower = 1> C;                    // number of cycles
  int<lower = 1> NN;                   // number of students in the current cycle
  matrix[N, M] x;                      // covariates
  vector[N] y;                         // outcomes
  int<lower = 1, upper = C> cycle[N];  // cycles
  vector[M] mu;                        // prior mean of beta
  vector<lower = 0>[M] sigma;          // prior sd of beta
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;
  
  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;
  
  vector[M] mu_std;
  vector<lower = 0>[M] sigma_std;
  
  // x[, 1] is the intercept
  x_std[, 1] = x[, 1];
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
  
  mu_std[1] = mu[1] - mu_y;
  sigma_std[1] = sigma[1] ^ 2;
  for (m in 2:M) {
    real r = sd_x[m] / sd_y;
    mu_std[m] = r * mu[m];
    sigma_std[m] = r * sigma[m];
    mu_std[1] += mu_x[m] * mu[m];
    sigma_std[1] += (mu_x[m] * sigma[m]) ^ 2;
  }
  mu_std[1] /= sd_y;
  sigma_std[1] = sqrt(sigma_std[1]) / sd_y;
}

parameters {
  vector[M] beta_std;
  real<lower = 0> sigma_y_std;
}

model {
  beta_std ~ normal(mu_std, sigma_std);
  sigma_y_std ~ cauchy(0, 2);
  
  y_std ~ normal(x_std * beta_std, sigma_y_std);
}

generated quantities {
  vector[M] beta;
  real<lower = 0> sigma_y = sigma_y_std * sd_y;
  vector[NN] log_lik;
  
  beta[1] = sd_y * beta_std[1] + mu_y;
  for (m in 2:M) {
    beta[m] = sd_y / sd_x[m] * beta_std[m];
    beta[1] -= beta[m] * mu_x[m];
  }
  
  for (n in 1:NN)
    log_lik[n] = normal_lpdf(y[n] | x[n, ] * beta, sigma_y);
}
'

##### inf
my.blr.inf <- function(dat_combined, n.chains = 4, prior.mean = mu, prior.sd = sigma) {
  writeLines(modelstring.inf, con = 'modelINF.stan')
  data <- c(data.stan(subset(dat_combined, cycle == 1)), list(mu = prior.mean, sigma = prior.sd))
  stan('modelINF.stan', data = data, iter = 20000, chains = 4)
}

# power priors
modelstring.pp <- '
data {
  int<lower = 1> N;                    // number of students
  int<lower = 1> M;                    // number of covariates
  int<lower = 1> C;                    // number of cycles
  int<lower = 1> NN;                   // number of students in the current cycle
  matrix[N, M] x;                      // covariates
  vector[N] y;                         // outcomes
  int<lower = 1, upper = C> cycle[N];  // cycles
  vector<lower = 0, upper = 1>[C] a;   // weight for each cycle
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;
  
  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;
  
  // x[, 1] is the intercept
  x_std[, 1] = x[, 1];
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
}

parameters {
  vector[M] beta_std;
  real<lower = 0> sigma_y_std;
}

model {
  beta_std ~ normal(0, 100);
  sigma_y_std ~ cauchy(0, 2);
  
  for (n in 1:N)
    target += a[cycle[n]] * normal_lpdf(y_std[n] | x_std[n, ] * beta_std, sigma_y_std);
}

generated quantities {
  vector[M] beta;
  real<lower = 0> sigma_y = sigma_y_std * sd_y;
  vector[NN] log_lik;
  
  beta[1] = sd_y * beta_std[1] + mu_y;
  for (m in 2:M) {
    beta[m] = sd_y / sd_x[m] * beta_std[m];
    beta[1] -= beta[m] * mu_x[m];
  }
  
  for (n in 1:NN)
    log_lik[n] = normal_lpdf(y[n] | x[n, ] * beta, sigma_y);
}
'

my.pp <- function(dat_combined, n.chains = 4, a = 0.5) {
  writeLines(modelstring.pp, con = 'modelPP.stan')
  data <- data.stan(dat_combined)
  data$a <- c(1, rep(a, data$C - 1))
  stan('modelPP.stan', data = data, iter = 20000, chains = 4)
}

# BDB
modelstring.bdb <- '
data {
  int<lower = 1> N;                    // number of students
  int<lower = 1> M;                    // number of covariates
  int<lower = 1> C;                    // number of cycles
  int<lower = 1> NN;                   // number of students in the current cycle
  matrix[N, M] x;                      // covariates
  vector[N] y;                         // outcomes
  int<lower = 1, upper = C> cycle[N];  // cycles
  real<lower = 0> nu1;                 // prior: tau2_beta ~ inv_gamma(nu1, nu2)
  real<lower = 0> nu2;                 // prior: tau2_beta ~ inv_gamma(nu1, nu2)
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;
  
  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;
  
  // x[, 1] is the intercept
  x_std[, 1] = x[, 1];
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
}

parameters {
  vector[M] beta_std[C];
  vector[M] mu_std;
  vector<lower = 0>[M] tau2_beta;
  
  real<lower = 0> sigma_y_std;
}

model {
  for (c in 1:C)
    beta_std[c] ~ normal(mu_std, sqrt(tau2_beta));
  mu_std ~ normal(0, 100);
  tau2_beta ~ inv_gamma(nu1, nu2);
  
  for (n in 1:N)
    y_std[n] ~ normal(x_std[n, ] * beta_std[cycle[n]], sigma_y_std);
  sigma_y_std ~ cauchy(0, 2);
}

generated quantities {
  vector[M] beta;
  real<lower = 0> sigma_y = sigma_y_std * sd_y;
  vector[NN] log_lik;
  
  beta[1] = sd_y * beta_std[1, 1] + mu_y;
  for (m in 2:M) {
    beta[m] = sd_y / sd_x[m] * beta_std[1, m];
    beta[1] -= beta[m] * mu_x[m];
  }
  
  for (n in 1:NN)
    log_lik[n] = normal_lpdf(y[n] | x[n, ] * beta, sigma_y);
}
'

my.bdb <- function(dat_combined, n.chains = 4, nu1 = 0.001, nu2 = 0.001) {
  writeLines(modelstring.bdb, con = 'modelBDB.stan')
  data <- data.stan(dat_combined)
  stan('modelBDB.stan', data = data, iter = 20000, chains = 4)
}



extract.stan <- function(fit, M = 5, coef = cf) {
  est <- summary(fit)$summary
  beta <-  est[grep('^beta\\[', rownames(est)), ][1:M, ]
  list(beta = beta,
       waic = waic(extract_log_lik(fit))$estimates[2:3, 1],
       looic = loo(fit)$estimates[2:3, 1],
       MSE = MSE(coef, beta[, 1]),
       bias = bias(coef, beta[, 1]),
       abs.bias = abs.bias(coef, beta[, 1]),
       per.bias = per.bias(coef, beta[, 1]),
       abs.per.bias = abs.per.bias(coef, beta[, 1]))
}



extract.bdb <- function(fit, M = 5, coef = cf) {
  est <- summary(fit)$summary
  beta <- est[grep('^beta\\[', rownames(est)), ][1:M, ]
  list(beta = beta,
       mu_std = est[grep('^mu_std\\[', rownames(est)), 1],
       sigma_beta_std = matrix(est[grep('^sigma_beta\\[', rownames(est)), 1], nrow = M),
       waic = waic(extract_log_lik(fit))$estimates[2:3, 1],
       looic = loo(fit)$estimates[2:3, 1],
       MSE = MSE(coef, beta[, 1]),
       bias = bias(coef, beta[, 1]),
       abs.bias = abs.bias(coef, beta[, 1]),
       per.bias = per.bias(coef, beta[, 1]),
       abs.per.bias = abs.per.bias(coef, beta[, 1]))
}
