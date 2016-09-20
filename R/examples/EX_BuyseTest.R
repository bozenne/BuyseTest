# reset the default value of the number of bootstrap samples
BuyseTest.options(n.bootstrap = 0) # no bootstrap

#### simulate some data ####
df.data <- simulBT(1e2, n.strata = 2)

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
                n.bootstrap = 1e3)
}
\dontshow{
  BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
                  n.bootstrap = 1e1, trace = 0)
}
summary(BT)
summary(BT, statistic = "winRatio") # win Ratio

## parallel boostrap
\dontrun{
  BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
                  n.bootstrap = 1e3, cpus = 2)
}
\dontshow{
  BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
                  n.bootstrap = 1e1, cpus = 2, trace = 0)
}
summary(BT)

## method Gehan is much faster but does not optimally handle censored observations
\dontrun{
  BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
                  n.bootstrap = 1e3, method = "Gehan")
}
\dontshow{
  BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
                  n.bootstrap = 1e1, method = "Gehan", trace = 0)
}
summary(BT)

#### one time to event endpoint: only differences in survival over 1 count ####
BT <- BuyseTest(Treatment ~ TTE(eventtime, threshold = 1, censoring = status), data=df.data)
summary(BT)

#### one time to event endpoint with a strata variable
BT <- BuyseTest(Treatment ~ strata + TTE(eventtime, censoring = status), data=df.data)
summary(BT)

#### several endpoints with a strata variable
BT <- BuyseTest(Treatment ~ strata + TTE(eventtime, 1, status) + Bin(toxicity) + TTE(eventtime, 0.5, status) + Cont(score, 1) + TTE(eventtime, 0.25, status), 
                data=df.data)
summary(BT)

#### real example : Veteran dataset of the survival package ####
#### Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#### method = "Gehan"

if(require(survival)){
\dontrun{
  data(veteran,package="survival")
  n.bootstrap <- 100
 
  ## method = "Gehan"
  BT_Gehan <- BuyseTest(trt ~ celltype + TTE(time,threshold=0,censoring=status), 
                        data=veteran, n.bootstrap=n.bootstrap,method="Gehan")
  
  summary_Gehan <- summary(BT_Gehan)
  summary_Gehan <- summary(BT_Gehan, statistic = "winRatio")
  
  ## method = "Peron"
  BT_Peron <- BuyseTest(data=veteran,endpoint="time",treatment="trt",strata="celltype",
                        type="timeToEvent",censoring="status",threshold=0,
                        n.bootstrap=n.bootstrap,method="Peron")
  
  summary_Peron <- summary(BT_Peron)
}
}
