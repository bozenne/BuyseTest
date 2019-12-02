### causal.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 22 2019 (09:53) 
## Version: 
## Last-Updated: nov 22 2019 (09:53) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## ** Causal
library(data.table)
library(BuyseTest)
library(pbapply)
library(lava)

m <- lvm(Y[1:1] ~ a*tau,
         X[0:1] ~ a*tau)
latent(m) <- ~tau

warper <- function(n, a, i){
    d <- sim(m, n = n, p = c("a" = a), latent = FALSE)
    idY <- sample.int(n/2)
    dY <- data.frame(score = d[idY,"Y"], group = "Y")
    dX <- data.frame(score = d[-idY,"X"], group = "X")
    dBT <- rbind(dX,dY)

    e.lm <- lm(score ~ group, data = dBT)
    e.BT <- BuyseTest(group ~ cont(score), data = dBT,
                      method.inference = "u-statistic",
                      trace = 0)
    row.ate <- data.frame(truth = mean(d$Y - d$X),
                          estimate = summary(e.lm)$coef["groupY","Estimate"],
                          se = summary(e.lm)$coef["groupY","Std. Error"],
                          lower.ci = confint(e.lm)["groupY","2.5 %"],
                          upper.ci = confint(e.lm)["groupY","97.5 %"],
                          p.value = summary(e.lm)$coef["groupY","Pr(>|t|)"],
                          statistic = "ate")
    GY0 <- (pnorm(d$X, mean = 0, sd = sqrt(1 + a^2))+pnorm(d$X, mean = 1, sd = sqrt(1 + a^2)))/2
    GY1 <- (pnorm(d$Y, mean = 0, sd = sqrt(1 + a^2))+pnorm(d$Y, mean = 1, sd = sqrt(1 + a^2)))/2
    row.quantile <- data.frame(truth = 2*mean(GY1 - GY0),
                          estimate = NA,
                          se = NA,
                          lower.ci = NA,
                          upper.ci = NA,
                          p.value = NA,
                          statistic = "quantile")
    row.benefit <- data.frame(truth = 2*mean(d$Y > d$X) - 1,
                              confint(e.BT),
                              statistic = "netBenefit")

    return(cbind(rbind(row.ate,
                       row.quantile,
                       row.benefit), iteration = i, cor = cor(d$X,d$Y)))
}

ls.sim <- pblapply(1:100, function(i){warper(n=500, a = 0, i = i)})
dt.sim <- as.data.table(do.call(rbind, ls.sim))
dt.sim[, .(truth = mean(truth), estimate = mean(estimate), cor = mean(cor)), by = "statistic"]


######################################################################
### causal.R ends here
