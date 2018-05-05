# reset the default value of the number of permuation sample
BuyseTest.options(method.inference = "none") # no permutation test

#### simulate some data ####
df.data <- simBuyseTest(1e2, n.strata = 2)

                                        # display 
if(require(survival)){
    resKM_tempo <- survfit(Surv(eventtime,status)~Treatment, data = df.data)
    plot(resKM_tempo)
}

#### one time to event endpoint ####
BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data)

summary(BT) # net chance in favor of treatment
summary(BT, percentage = FALSE)  
summary(BT, statistic = "winRatio") # win Ratio

## bootstrap to compute the CI
\dontrun{
    BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
                    method.inference = "permutation", n.resampling = 1e3)
}
\dontshow{
    BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
                    method.inference = "permutation", n.resampling = 1e1, trace = 0)
}
summary(BT, statistic = "netChance") ## default
summary(BT, statistic = "winRatio") 

## parallel boostrap
\dontrun{
    BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
                    method.inference = "permutation", n.resampling = 1e3, cpus = 2)
    summary(BT)
}

## method Gehan is much faster but does not optimally handle censored observations
BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
                method.tte = "Gehan", trace = 0)
summary(BT)

#### one time to event endpoint: only differences in survival over 1 unit ####
BT <- BuyseTest(Treatment ~ TTE(eventtime, threshold = 1, censoring = status), data=df.data)
summary(BT)

#### one time to event endpoint with a strata variable
BT <- BuyseTest(Treatment ~ strata + TTE(eventtime, censoring = status), data=df.data)
summary(BT)

#### several endpoints with a strata variable
f <- Treatment ~ strata + T(eventtime, 1, status) + B(toxicity) 
f <- update(f, 
            ~. + T(eventtime, 0.5, status) + C(score, 1) + T(eventtime, 0.25, status))

BT <- BuyseTest(f, data=df.data)
summary(BT)

#### real example : Veteran dataset of the survival package ####
#### Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#### method.tte = "Gehan"

if(require(survival)){
\dontrun{
  data(veteran,package="survival")
 
  ## method.tte = "Gehan"
  BT_Gehan <- BuyseTest(trt ~ celltype + TTE(time,threshold=0,censoring=status), 
                        data=veteran, method.tte="Gehan",
                        method.inference = "permutation", n.resampling = 1e3)
  
  summary_Gehan <- summary(BT_Gehan)
  summary_Gehan <- summary(BT_Gehan, statistic = "winRatio")
  
  ## method.tte = "Peron"
  BT_Peron <- BuyseTest(trt ~ celltype + TTE(time,threshold=0,censoring=status), 
                        data=veteran, method.tte="Peron",
                        method.inference = "permutation", n.resampling = 1e3)

  class(BT_Peron)
  summary(BT_Peron)
}
}
