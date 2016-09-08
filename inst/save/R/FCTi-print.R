#' @name internal-print
#' @rdname internal-print
#' @title internal functions for BuyseTest - display
#' @aliases printGeneral printBoostrap
#' 
#' @description Functions called by \code{\link{BuyseTest}} to display what is going on.
#' 
#' @details
#' \code{printGeneral} display general settings. \cr
#' \code{printBoostrap} display bootstrap settings. \cr
#' 
#' 
#' @keywords function internal BuyseTest


printGeneral <- function(levels.treatment,
                         levels.strata, n.strata,
                         endpoint, threshold, censoring, type, D, D.TTE,
                         method, Wscheme, threshold_TTEM1){
  cat("Settings (general) \n")
  cat("   # chosen reference : Control = ",levels.treatment[1]," and Treatment = ",levels.treatment[2],"\n")
  cat("   # number of endpoints : ",D," \n")
  iter_mTTE <- 0
  for (iter_m in 1:D) {
    cat("      >  endpoint ", iter_m, " = \"", endpoint[iter_m], "\" - type = \"", c("binary", "continuous", "timeToEvent")[type[iter_m]], "\"", sep = "")
    if (type[iter_m] %in% c(2,3)) {cat(" | threshold =",threshold[iter_m])}
    if (type[iter_m] == 3) {iter_mTTE <- iter_mTTE + 1 ; cat(" | censoring = \"", censoring[iter_mTTE], "\"", sep = "");}
    cat("\n")
  }
  cat("   # n.strata = ", n.strata, " : ", paste(levels.strata,collapse = " "), "\n")
  cat("   # management of censored survival pairs : ")
  switch(method,
         "Gehan" = cat("uninformative pairs \n"),
         "Peto" = cat("imputation using one survival curve estimated on all patients \n"),
         "Efron" = cat("imputation using different survival curve for control and treatment patients \n"),
         "Peron" = cat("imputation using different survival curve for control and treatment patients \n")
  ) 
  if (method %in% c("Peto","Efron","Peron")) {
    
    cat("   # weights of the pairs relatively to the enpoints : \n")
    print(Wscheme)
    
    cat("   # intervals thresholds for survival endpoints : \n")
    
    if (length(threshold_TTEM1) > 0) {
      threshold_TTEM1.display <- threshold_TTEM1
      threshold_TTEM1.display[threshold_TTEM1.display < 0] <- +Inf
    }else{
      threshold_TTEM1.display <- +Inf
    }
    
    threshold.display <- rbind(sapply(1:D.TTE,
                                      function(x){paste(c("[",round(threshold_TTEM1.display[x],4),
                                                          " ; ",round(threshold[type == 3][x],4),
                                                          "] "), collapse = "")}))
    colnames(threshold.display) <- endpoint[type == 3]      
    rownames(threshold.display) <- "threshold interval"
    print(threshold.display)
  }
}

printBoostrap <- function(prob.alloc, n.bootstrap, stratified, cpus, time, seed){
  
  if (time[3] == 0) {
    time.punctual <- "<0.001 s"
    time.bootstrap <- paste0("<",signif(0.001*n.bootstrap/cpus,4)," s")
  }else{
    time.punctual <- paste(time[3],"s")
    time.bootstrap <- paste(signif(time[3]*n.bootstrap/cpus,4),"s")
  }
  
  cat("Settings (", if (stratified) {"stratified "}, "bootstrap) \n",
      "   # resampling probability for assignement to the treatment group: ", prob.alloc, "\n",
      "   # requested time for one sample: ", time.punctual, "\n",
      "   # estimated time for ", n.bootstrap, " samples with ", cpus, " core", if (cpus > 1) {"s"}, ": ", time.bootstrap, "\n", sep = "")
  if (!is.null(seed)) {
    cat("   # seed", if (cpus > 1) {"s"}, ": ",paste(seq(seed,seed + cpus - 1), collapse = " "), " \n", sep = "")       
  }
}