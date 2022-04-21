#Single level case study codes
source("Functions.R")
all.dt <- read.csv("multilevel.all.cycles.csv")
all.dt[,c(1:5,7)] <- NULL
head(all.dt)

get.beta.loo <- function(fit, M = 5) {
  est <- summary(fit)$summary
  beta <- est[grep('^beta\\[', rownames(est)), ][1:M, ]
  list(beta = beta,
       looic = loo(fit))
}

#prepare all.dt
blr.non.inf <- my.blr(all.dt)
blr.non.inf.res <- get.beta.loo(blr.non.inf)

#inf
pri.mn <- c(407.385619, -11.418260 ,  5.149994,  25.584788,   3.617521)
pri.sd <- c(6.652371, 2.253203, 0.493518, 1.189917, 3.069737 )

blr.inf <- my.blr.inf(all.dt, prior.mean = pri.mn, prior.sd = pri.sd)
blr.inf.res <- get.beta.loo(blr.inf)

#full.borrowing - pooling
blr.pool <- my.blr.pool(all.dt)
blr.pool.res <- get.beta.loo(blr.pool)

#BDB
bdb_.001 <- my.bdb(all.dt, nu1 = .001, nu2 = .001)
bdb_.001.res <- get.beta.loo(bdb_.001)

bdb_.01 <- my.bdb(all.dt, nu1 = .01, nu2 = .01)
bdb_.01.res <- get.beta.loo(bdb_.01)

bdb_.1 <- my.bdb(all.dt, nu1 = .1, nu2 = .1)
bdb_.1.res <- get.beta.loo(bdb_.1)

bdb_1_1 <- my.bdb(all.dt, nu1 = 1, nu2 = 1)
bdb_1_1.res <- get.beta.loo(bdb_1_1)

bdb_1_.1 <- my.bdb(all.dt, nu1 = 1, nu2 = .1)
bdb_1_.1.res <- get.beta.loo(bdb_1_.1)

bdb_1_.01 <- my.bdb(all.dt, nu1 = 1, nu2 = .01)
bdb_1_.01.res <- get.beta.loo(bdb_1_.01)

bdb_1_.001 <- my.bdb(all.dt, nu1 = 1, nu2 = .001)
bdb_1_.001.res <- get.beta.loo(bdb_1_.001)

pp.25 <- my.pp(all.dt, a = 0.25)
pp.25.res <- get.beta.loo(pp.25)

pp.50 <- my.pp(all.dt, a = 0.5)
pp.50.res <- get.beta.loo(pp.50)

pp.75 <- my.pp(all.dt, a = 0.75)
pp.75.res <- get.beta.loo(pp.75)

save.image("case_study.RData")

est <- rbind(blr.non.inf.res$beta,
             blr.inf.res$beta,
             blr.pool.res$beta,
             bdb_.001.res$beta,
             bdb_.01.res$beta,
             bdb_.1.res$beta,
             bdb_1_1.res$beta,
             bdb_1_.1.res$beta,
             bdb_1_.01.res$beta,
             bdb_1_.001.res$beta,
             bdb_1_3.res$beta,
             pp.0.res$beta,
             pp.25.res$beta,
             pp.50.res$beta,
             pp.75.res$beta,
             pp.1.res$beta)

blr.non.inf.res$looic$se_looic

loo <- rbind(blr.non.inf.res$looic,
             blr.inf.res$looic,
             blr.pool.res$looic,
             bdb_.001.res$looic,
             bdb_.01.res$looic,
             bdb_.1.res$looic,
             bdb_1_1.res$looic,
             bdb_1_.1.res$looic,
             bdb_1_.01.res$looic,
             bdb_1_.001.res$looic,
             bdb_1_3.res$looic,
             pp.0.res$looic,
             pp.25.res$looic,
             pp.50.res$looic,
             pp.75.res$looic,
             pp.1.res$looic)

loo_se <- rbind(blr.non.inf.res$looic$se_looic,
                blr.inf.res$looic$se_looic,
                blr.pool.res$looic$se_looic,
                bdb_.001.res$looic$se_looic,
                bdb_.01.res$looic$se_looic,
                bdb_.1.res$looic$se_looic,
                bdb_1_1.res$looic$se_looic,
                bdb_1_.1.res$looic$se_looic,
                bdb_1_.01.res$looic$se_looic,
                bdb_1_.001.res$looic$se_looic,
                pp.25.res$looic$se_looic,
                pp.50.res$looic$se_looic,
                pp.75.res$looic$se_looic)

row.names(loo) <- c("blr.non.inf.res", "blr.inf.res", "blr.pool.res", "bdb_.001.res", "bdb_.01.res", 
                    "bdb_.1.res", "bdb_1_1.res", "bdb_1_.1.res", "bdb_1_.01.res", "bdb_1_.001.res", "bdb_1_3.res", "pp.0.res", 
                    "pp.25.res", "pp.50.res", "pp.75.res", "pp.1.res")

library(xlsx)

write.xlsx(est, paste0("CaseResults.xlsx"), sheetName = "Estimates", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

write.xlsx(loo, paste0("CaseResults.xlsx"), sheetName = "Loo", 
           col.names = TRUE, row.names = TRUE, append = TRUE)

