## * Single Wilcoxon test

## chunk 2
x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
df <- rbind(data.frame(value = x, group="x"),
            data.frame(value = y, group="y"))

## chunk 3
wilcox.test(value ~ group, data = df)

## chunk 4
library(asht)
wmwTest(value ~ group, data = df, method = "asymptotic")

## chunk 5
wmwTest(value ~ group, data = df, method = "exact.ce")

## chunk 6
eperm.BT <- BuyseTest(group ~ cont(value), data = df, add.halfNeutral = TRUE,
                      method.inference = "permutation", n.resampling = 10000,
                      trace = FALSE, seed = 10)
confint(eperm.BT, statistic = "favorable")

## chunk 7
BuyseTest.options(order.Hprojection=2)
eU.BT <- BuyseTest(group ~ cont(value), data = df,
                  method.inference = "u-statistic",
                  add.halfNeutral = TRUE, trace = FALSE)
confint(eU.BT, statistic = "favorable")

## chunk 8
etperm.BT <- BuyseTest(group ~ cont(value), data = df, add.halfNeutral = TRUE,
                      method.inference = "studentized permutation", n.resampling = 10000,
                      trace = FALSE, seed = 10)
confint(etperm.BT, statistic = "favorable")

## * Multiple Wilcoxon tests

## chunk 9
set.seed(35)
dt <- simBuyseTest(n.T=25, n.strata = 5)
dt$id <- paste0("id",1:NROW(dt))
dt$strata <- as.character(dt$strata) 
head(dt)

## chunk 10
BuyseTest.options(order.Hprojection=1);BuyseTest.options(trace=0)

ls.BT <- list("b-a=0" = BuyseTest(strata ~ cont(score), add.halfNeutral = TRUE,
                                  data = dt[dt$strata %in% c("a","b"),],
                                  method.inference = "u-statistic"),
              "c-a=0" = BuyseTest(strata ~ cont(score), add.halfNeutral = TRUE,
                                  data = dt[dt$strata %in% c("a","c"),],
                                  method.inference = "u-statistic"),
              "d-a=0" = BuyseTest(strata ~ cont(score), add.halfNeutral = TRUE,
                                  data = dt[dt$strata %in% c("a","d"),],
                                  method.inference = "u-statistic"),
              "e-a=0" = BuyseTest(strata ~ cont(score), add.halfNeutral = TRUE,
                                  data = dt[dt$strata %in% c("a","e"),],
                                  method.inference = "u-statistic")
              )

M.confint <- do.call(rbind,lapply(ls.BT,confint, statistic = "favorable"))
cbind(M.confint,adj.p.value = p.adjust(M.confint[,"p.value"], method = "bonferroni"))

## chunk 11
e.mc <- BuyseMultComp(ls.BT, statistic = "favorable", cluster = "id", global = TRUE)
print(e.mc, cols = c("estimate","se","p.value","adj.p.value"))

## chunk 12
M.cor <- cov2cor(crossprod(e.mc$iid))
dimnames(M.cor) <- list(names(ls.BT),names(ls.BT))
M.cor

