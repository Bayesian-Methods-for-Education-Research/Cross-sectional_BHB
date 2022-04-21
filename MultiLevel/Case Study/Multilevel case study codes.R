#Multilevel case study codes
library(tidyverse)
library(rstan)
library(loo)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
parallel:::setDefaultClusterOptions(setup_strategy = 'sequential')

set.seed(98765)
set.seed(round(runif(1, 0, 1e8)))

# combine all csv data files
all <- data.frame()
for (file in sprintf('ps%02d.csv', seq(3, 18, by = 3), '.csv'))
  all <- rbind(all, read.csv(file))
C <- max(all$cycle)

# drop schools with less than 10 students
n.sch <- group_by(all, cycle, SCHOOLID) %>%
  summarize(n = n())
data <- inner_join(all, n.sch, by = c('cycle', 'SCHOOLID')) %>%
  filter(n >= 10) %>%
  select(-n) %>%
  rename(schid = SCHOOLID)

# outcome model: please change it if needed!
# MM is the number of level 1 parameters that have random effects (including the intercept)
# make sure that variables 1:MM are those have random effects (that is, put those variables first)
# case 1: allow all interactions
model <- ~ Female + PARED + HOMEPOS + IMMIG + TCSHORT + STRATIO + Female:TCSHORT
MM <- 2
# case 2: only the slope of HOMEPOS is explained by level 2 variables
# model <- ~ HOMEPOS * (TCSHORT + STRATIO) + Female + PARED + IMMIG
# MM <- 2

# fit noninformative prior models to historical cycles to find the informative prior for the current cycle
estimate <- list()
fits <- list()
rhat <- rep(NA, C)
for (c in 2:C) {
  print(c)
  dat <- filter(data, cycle == c) %>%
    mutate(schid = as.integer(as.factor(schid)), cycle = 1)
  x <- as.data.frame(model.matrix(model, dat))
  fit <- stan('reg.stan', chains = 4, iter = 30000, thin = 10, save_warmup = F,
              data = list(N = nrow(dat), NN = 1, M = ncol(x), MM = MM, S = max(dat$schid), C = 1, x = x, y = dat$PV1MATH, sch = dat$schid, cycle = dat$cycle))
  est <- summary(fit)$summary
  rhat[c] <- max(est[, 'Rhat'], na.rm = T)
  name <- rownames(est)
  beta <- est[grep('^beta\\[', name), c(1, 3)]
  names(beta) <- colnames(x)[1:nrow(beta)]
  print(fit, 'beta')
  estimate[[c]] <- list(beta = beta[, 1], sigma.beta = beta[, 2],
                        sigma.e = est[grep('^sigma_e$', name), 1],
                        #tau.u = est[grep('^tau_u\\[', name)],
                        sigma.u = as.matrix(est[grep('^sigma_u\\[', name), 1], nrow = 1, byrow = T))
  fits[[c]] <- fit
}
# check convergence
print(rhat)

coefs <- do.call(rbind, lapply(estimate, unlist))
betas <- coefs[, 1:ncol(x)]
sigmas.beta <- coefs[, (ncol(x) + 1):(ncol(x) * 2)]
sigmas <- coefs[, (ncol(x) * 2 + 1):ncol(coefs)]
save(list = c('data', 'fits', 'betas', 'sigmas.beta', 'sigmas'), file = paste0('input.RData'))

load('input.RData')

methods <- c('reg_noninf', 'reg_pooling', 'inf', 'pp_0.25', 'pp_0.5', 'pp_0.75',
             sort(unlist(lapply(c(1, 20), function(nu) {
               paste0(c('bdb_0.001_0.001', 'bdb_0.01_0.01', 'bdb_0.1_0.1', 'bdb_1_1', 'bdb_1_0.001', 'bdb_1_0.01', 'bdb_1_0.1'),
                      '_', nu)
             }))))
all <- data

output <- list()
for (method in 1:length(methods)) {
  print(methods[method])
  s <- strsplit(methods[method], '_')[[1]]
  data <- if (s[2] == 'noninf') subset(all, cycle == 1) else all[rank(all$cycle, ties.method = 'first'), ]
  data$schid <- as.integer(as.factor(data$schid))
  x <- as.data.frame(model.matrix(model, data))
  dat <- list(N = nrow(x), NN = sum(data$cycle == 1), M = ncol(x), MM = MM,
              S = max(data$schid), C = max(data$cycle), x = x, y = as.numeric(data$PV1MATH),
              sch = data$schid, cycle = data$cycle, cycle_sch = unlist(by(data$cycle, data$schid, median, simplify = F)))
  if (s[1] == 'inf') {
    dat$mu <- colMeans(betas[-1, ])
    dat$sigma <- colMeans(sigmas.beta[-1, ])
  } else if (s[1] == 'pp') {
    dat$a <- c(1, rep(as.numeric(s[2]), dat$C - 1))
  } else if (s[1] == 'bdb') {
    dat$nu1 <- as.numeric(s[2])
    dat$nu2 <- as.numeric(s[3])
    dat$nu <- as.numeric(s[4])
  }            
  
  fit <- stan(paste0(s[1], '.stan'), data = dat, iter = 30000, chains = 4, thin = 10,
              save_warmup = F, pars = c('beta', 'sigma_e', 'sigma_u', 'log_lik'), include = T)
  est <- summary(fit)$summary
  ic <- list(waic = waic(extract_log_lik(fit))$estimates, looic = loo(fit)$estimates)
  max.rhat <- max(est[, 'Rhat'], na.rm = T)
  if (max.rhat > 1.1)
    print(paste0('Error: ', method, ', ', max.rhat))
  print(fit, 'beta')
  print(max.rhat)
  output[[method]] <- list(fit = fit, est = est, ic = ic, method = methods[method])
  rm(list = c('fit', 'est', 'ic'))
  gc(F)
}

save(list = c('output'), file = 'output.RData')
