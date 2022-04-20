#Single level simulation study codes

a <- "cond1" #can go up to cond10
source("Functions.R")

rep=1 #1:500
set.seed(2926100)
seed.number <- round(runif(1000,10000000, 99999999),0)
my.seed <- seed.number[rep]
set.seed(my.seed)

cond <- read.csv("cond.100.csv") #conditions    - read *n.500.cvs or *.n.2000.csv for other sample size conditions 
ps03 <- read.csv("ps03.n100.csv") #scaled data - read *n.500.cvs or *.n.2000.csv for other sample size conditions 
ps06 <- read.csv("ps06.n100.csv") #scaled data - read *n.500.cvs or *.n.2000.csv for other sample size conditions 
ps09 <- read.csv("ps09.n100.csv") #scaled data - read *n.500.cvs or *.n.2000.csv for other sample size conditions 
ps12 <- read.csv("ps12.n100.csv") #scaled data - read *n.500.cvs or *.n.2000.csv for other sample size conditions 
ps15 <- read.csv("ps15.n100.csv") #scaled data - read *n.500.cvs or *.n.2000.csv for other sample size conditions 
ps18 <- read.csv("ps18.csv") #not scaled
ps18[,c(1:5,7,12)] <- NULL
ps03[,c(1:5,7)] <- NULL
ps06[,c(1:5,7)] <- NULL
ps09[,c(1:5,7)] <- NULL
ps12[,c(1:5,7)] <- NULL
ps15[,c(1:5,7)] <- NULL

cf <- cond[,a]

ps18.syn <- generate.pisa(ps18, n =100, coef=cf, res.std.err = 81) #generate non-scaled - change n = 100 accordingly

#prepare all.dt
ps18.syn$cycle <- 1
all.dt <- rbind(ps03, ps06, ps09, ps12, ps15, ps18.syn)

#uninformatiove
blr.non.inf <- my.blr(all.dt)
est.blr.non.inf <- extract.stan(blr.non.inf)

#inf
pri.mn <- cond[,"cond5"]
pri.sd <- cond[,"cf.sd"]
blr.inf <- my.blr.inf(all.dt, prior.mean = pri.mn, prior.sd = pri.sd)
est.blr.inf <- extract.stan(blr.inf)

#full.borrowing - pooling
blr.pool <- my.blr.pool(all.dt)
est.blr.pool <- extract.stan(blr.pool)

#BDB
bdb_.001 <- my.bdb(all.dt, nu1 = .001, nu2 = .001)
bdb_.01 <- my.bdb(all.dt, nu1 = .01, nu2 = .01)
bdb_.1 <- my.bdb(all.dt, nu1 = .1, nu2 = .1)
bdb_1_1 <- my.bdb(all.dt, nu1 = 1, nu2 = 1)
bdb_1_.1 <- my.bdb(all.dt, nu1 = 1, nu2 = .1)
bdb_1_.01 <- my.bdb(all.dt, nu1 = 1, nu2 = .01)
bdb_1_.001 <- my.bdb(all.dt, nu1 = 1, nu2 = .001)
bdb_1_3 <- my.bdb(all.dt, nu1 = 1, nu2 = 3)

est.bdb_.001 <- extract.bdb(bdb_.001)
est.bdb_.01 <- extract.bdb(bdb_.01)
est.bdb_.1 <- extract.bdb(bdb_.1)
est.bdb_1_1 <- extract.bdb(bdb_1_1)
est.bdb_1_.1 <- extract.bdb(bdb_1_.1)
est.bdb_1_.01 <- extract.bdb(bdb_1_.01)
est.bdb_1_.001 <- extract.bdb(bdb_1_.001)
est.bdb_1_3 <- extract.bdb(bdb_1_3)

# pp
pp.0 <- my.pp(all.dt, a = 0)
pp.25 <- my.pp(all.dt, a = 0.25)
pp.50 <- my.pp(all.dt, a = 0.5)
pp.75 <- my.pp(all.dt, a = 0.75)
pp.1 <- my.pp(all.dt, a = 1)

est.pp.0 <- extract.stan(pp.0)
est.pp.25 <- extract.stan(pp.25)
est.pp.50 <- extract.stan(pp.50)
est.pp.75 <- extract.stan(pp.75)
est.pp.1 <- extract.stan(pp.1)

#here are the results

rmse <- cbind(est.blr.non.inf$MSE, est.blr.inf$MSE, est.blr.pool$MSE, est.bdb_.001$MSE, est.bdb_.01$MSE, 
              est.bdb_.1$MSE, est.bdb_1_1$MSE, est.bdb_1_.1$MSE, est.bdb_1_.01$MSE, est.bdb_1_.001$MSE, 
              est.bdb_1_3$MSE, est.pp.0$MSE, est.pp.25$MSE, est.pp.50$MSE, est.pp.75$MSE, est.pp.1$MSE)

rbias <- cbind(est.blr.non.inf$bias, est.blr.inf$bias, est.blr.pool$bias, est.bdb_.001$bias, est.bdb_.01$bias, 
               est.bdb_.1$bias, est.bdb_1_1$bias, est.bdb_1_.1$bias, est.bdb_1_.01$bias, est.bdb_1_.001$bias, 
               est.bdb_1_3$bias, est.pp.0$bias, est.pp.25$bias, est.pp.50$bias, est.pp.75$bias, est.pp.1$bias)

rabias <- cbind(est.blr.non.inf$abs.bias, est.blr.inf$abs.bias, est.blr.pool$abs.bias, est.bdb_.001$abs.bias, est.bdb_.01$abs.bias, 
                est.bdb_.1$abs.bias, est.bdb_1_1$abs.bias, est.bdb_1_.1$abs.bias, est.bdb_1_.01$abs.bias, est.bdb_1_.001$abs.bias, 
                est.bdb_1_3$abs.bias, est.pp.0$abs.bias, est.pp.25$abs.bias, est.pp.50$abs.bias, est.pp.75$abs.bias, est.pp.1$abs.bias)

rpbias <- cbind(est.blr.non.inf$per.bias, est.blr.inf$per.bias, est.blr.pool$per.bias, est.bdb_.001$per.bias, est.bdb_.01$per.bias, 
                est.bdb_.1$per.bias, est.bdb_1_1$per.bias, est.bdb_1_.1$per.bias, est.bdb_1_.01$per.bias, est.bdb_1_.001$per.bias, 
                est.bdb_1_3$per.bias, est.pp.0$per.bias, est.pp.25$per.bias, est.pp.50$per.bias, est.pp.75$per.bias, est.pp.1$per.bias)

rapbias <- cbind(est.blr.non.inf$abs.per.bias, est.blr.inf$abs.per.bias, est.blr.pool$abs.per.bias, est.bdb_.001$abs.per.bias, 
                 est.bdb_.01$abs.per.bias, est.bdb_.1$abs.per.bias, est.bdb_1_1$abs.per.bias, est.bdb_1_.1$abs.per.bias, 
                 est.bdb_1_.01$abs.per.bias, est.bdb_1_.001$abs.per.bias, est.bdb_1_3$abs.per.bias, est.pp.0$abs.per.bias, 
                 est.pp.25$abs.per.bias, est.pp.50$abs.per.bias, est.pp.75$abs.per.bias, est.pp.1$abs.per.bias)

rw <- cbind(est.blr.non.inf$waic, est.blr.inf$waic, est.blr.pool$waic, est.bdb_.001$waic, est.bdb_.01$waic, est.bdb_.1$waic,
            est.bdb_1_1$waic, est.bdb_1_.1$waic, est.bdb_1_.01$waic, est.bdb_1_.001$waic, est.bdb_1_3$waic, est.pp.0$waic, 
            est.pp.25$waic, est.pp.50$waic, est.pp.75$waic, est.pp.1$waic)

rl <- cbind(est.blr.non.inf$looic, est.blr.inf$looic, est.blr.pool$looic, est.bdb_.001$looic, est.bdb_.01$looic, est.bdb_.1$looic, 
            est.bdb_1_1$looic, est.bdb_1_.1$looic, est.bdb_1_.01$looic, est.bdb_1_.001$looic, est.bdb_1_3$looic, est.pp.0$looic, 
            est.pp.25$looic, est.pp.50$looic, est.pp.75$looic, est.pp.1$looic)

rest <- rbind(est.blr.non.inf$beta, est.blr.inf$beta, est.blr.pool$beta, est.bdb_.001$beta, est.bdb_.01$beta, est.bdb_.1$beta, 
              est.bdb_1_1$beta, est.bdb_1_.1$beta, est.bdb_1_.01$beta, est.bdb_1_.001$beta, est.bdb_1_3$beta, est.pp.0$beta, 
              est.pp.25$beta, est.pp.50$beta, est.pp.75$beta, est.pp.1$beta)

r1 <- as.vector(rmse)
r2 <- as.vector(rbias)
r3 <- as.vector(rabias)
r4 <- as.vector(rpbias)
r5 <- as.vector(rapbias)
r6 <- as.vector(rw)
r7 <- as.vector(rl)
r8 <- as.vector(rest)

res <- c(r1, r2, r3, r4, r5, r6, r7, r8)

write.table(t(res), paste0("results.100.",a,".",rep,".csv"), row.names = FALSE, col.names = FALSE)
