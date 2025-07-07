## chunk 2
data(diabetic, package = "survival")
head(diabetic)

## chunk 3
diabetic$juvenile <- diabetic$age <= 19
library(LMMstar)
summarize(age ~ juvenile, data = diabetic[!duplicated(diabetic$id),])

## chunk 4
diabeticJ <- diabetic[diabetic$juvenile,]

## * Wald methods (Gehan scoring rule)

## chunk 5
library(BuyseTest)
e.BTjuv <- BuyseTest(trt ~ tte(time,status) + strata(id, match = TRUE), 
                     data = diabeticJ, trace = FALSE,
                     scoring.rule =  "Gehan")
model.tables(e.BTjuv, percentage = FALSE)

## chunk 6
mean(coef(e.BTjuv, strata = TRUE))

## chunk 7
N <- nobs(e.BTjuv)["pairs"]
pw <- coef(e.BTjuv, statistic = "favorable")
pl <- coef(e.BTjuv, statistic = "unfavorable")
sqrt((pw + pl - (pw - pl)^2)/N)

## chunk 8
confint(e.BTjuv)

## chunk 9
confint(e.BTjuv, transform = FALSE)

## chunk 10
sqrt(var(coef(e.BTjuv, strata = TRUE))/N)

## chunk 11
sqrt(var(coef(e.BTjuv, strata = TRUE))/N) * sqrt((N-1)/N)

## * MOVER method (Gehan scoring rule)

## chunk 12
mover(e.BTjuv)

## * Wald methods (Peron scoring rule)

## chunk 13
library(prodlim)
e.BTjuv2 <- BuyseTest(trt ~ tte(time,status) + strata(id, match = TRUE), 
                      data = diabeticJ, trace = FALSE,
                      model.tte = prodlim(Hist(time,status)~ trt, data = diabeticJ))
model.tables(e.BTjuv2, percentage = FALSE)

## chunk 14
c(sqrt(var(coef(e.BTjuv2, strata = TRUE))/N),
  sqrt(var(coef(e.BTjuv2, strata = TRUE))/N) * sqrt((N-1)/N)
  )

## chunk 15
model.tte <- prodlim(Hist(time,status)~ trt, data = diabeticJ)
attr(model.tte, "iidNuisance") <- FALSE
confint(BuyseTest(trt ~ tte(time,status) + strata(id, match = TRUE), 
                  data = diabeticJ, trace = FALSE,
                  model.tte = model.tte))

## chunk 16
pw.peron <- coef(e.BTjuv2, statistic = "favorable")
pl.peron <- coef(e.BTjuv2, statistic = "unfavorable")
sqrt((pw.peron + pl.peron - (pw.peron - pl.peron)^2)/N)

## chunk 17
confint(e.BTjuv2)

## chunk 18
pw.peronS <- coef(e.BTjuv2, statistic = "favorable", strata = TRUE)
pl.peronS <- coef(e.BTjuv2, statistic = "unfavorable", strata = TRUE)
Hterm1 <- (pw.peronS - pl.peronS) - (pw.peron - pl.peron)

## chunk 19
Hterm2.obs <- e.BTjuv2@iidNuisance$favorable - e.BTjuv2@iidNuisance$unfavorable
Hterm2.pair <- Hterm2.obs[diabeticJ$trt==0]
table(Hterm2.obs[diabeticJ$trt==1])

## chunk 20
c(average = sqrt(crossprod((Hterm1/N))),
  nuisance = sqrt(crossprod((Hterm2.pair/N))),
  all = sqrt(crossprod((Hterm1/N + Hterm2.pair/N))))

## * References
