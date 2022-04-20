#Multilevel simulation study codes
library(rstan)
library(loo)
options(mc.cores = 4)

condition <- 1 #1 to 10
replication <- 1 #1 to 500

set.seed(condition ^ 2)
repeat {
  seeds <- unique(round(runif(3000, 0, 1e8)))
  if (length(seeds) >= 1000) break
}
set.seed(seeds[replication])

methods <- c('reg_noninf', 'reg_pooling', 'inf', 'pp_0.25', 'pp_0.5', 'pp_0.75',
             sort(unlist(lapply(c(1, 20), function(nu) {
               paste0(c('bdb_0.001_0.001', 'bdb_0.01_0.01', 'bdb_0.1_0.1', 'bdb_1_1', 'bdb_1_0.001', 'bdb_1_0.01', 'bdb_1_0.1'),
                      '_', nu)
             }))))
load('input-multilevelSimulation_10sch-40Students.RData')

beta <- colMeans(betas[-1, ])
sigma <- colMeans(sigmas[-1, ])
beta <- c(beta[1],
          if (condition == 10)
            -beta[-1]
          else
            beta[-1] + abs(beta[-1]) * c(-0.8, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 0.8)[condition]
)

dat <- subset(data, cycle == 1)
u <- rnorm(max(dat$schid), 0, sqrt(sigma[2]))
dat$PV1MATH <- as.matrix(cbind(1, dat[1:(ncol(dat) - 3)])) %*% beta + u[dat$schid] + rnorm(nrow(dat), 0, sigma[1])

all <- rbind(dat, subset(data, cycle != 1))

output <- list()
for (method in 1:length(methods)) {
  print(methods[method])
  s <- strsplit(methods[method], '_')[[1]]
  data <- if (s[1] == 'inf' || s[2] == 'noninf')
    subset(all, cycle == 1)
  else
    all[rank(all$cycle, ties.method = 'first'), ]
  x <- as.data.frame(model.matrix(~ Female + PARED + HOMEPOS + IMMIG + TCSHORT + STRATIO, data))
  dat <- list(N = nrow(x), NN = sum(data$cycle == 1), M = ncol(x), MM = 1,
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
    write.table(c(method, max.rhat), paste0('error_', condition, '_', replication, '.log'), T)
  print(fit, 'beta')
  print(max.rhat)
  output[[method]] <- list(est = est, ic = ic, method = methods[method])
  rm(list = c('fit', 'est', 'ic'))
  gc(F)
}

save(list = c('beta', 'sigma', 'output'), file = paste0('output_', condition, '_', replication, '.RData'))


Additional R codes 