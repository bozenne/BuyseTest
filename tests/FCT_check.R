validPairs <- function(BuyseRes, type = c("strata","sum")){
  
  BuyseSummary <- summary(BuyseRes, show = NULL)
  enpoint_threshold <- paste(BuyseSummary$nb$endpoint,BuyseSummary$nb$threshold, sep = "_")
  endpoints <- unique(enpoint_threshold)
  D <- length(endpoints)
  index_strata <- which(BuyseSummary$nb$strata!="global")
  levels.strata <- unique(BuyseSummary$nb$strata[index_strata])
  n.strata <- length(levels.strata) 
  
  diff <- NULL
  
  if("strata" %in% type){
    diff.strata <- matrix(NA,nrow=(1+n.strata)*D,ncol=5)
    colnames(diff.strata) <- c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")
    rownames(diff.strata) <- paste(BuyseSummary$nb$endpoint, " th=",BuyseSummary$nbthreshold ," strata=", BuyseSummary$nb$strata)
    for(iter_endpoint in 1:D){
      index_endpoint <- which(enpoint_threshold==endpoints[iter_endpoint])  
      res_strata <- as.matrix(BuyseSummary$nb[intersect(index_strata,index_endpoint),c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")])
      diff.strata[intersect(index_strata,index_endpoint),] <- t(apply(res_strata,1,function(x){x-res_strata[1,,drop=TRUE]}))
    }
    diff <- cbind(diff, diff.strata)
  }
  if("sum" %in% type){
    diff.sum <- cbind(sum = BuyseSummary$nb$n.total - (BuyseSummary$nb$n.favorable+BuyseSummary$nb$n.unfavorable+BuyseSummary$nb$n.neutral+BuyseSummary$nb$n.uninf))
    rownames(diff.sum) <- paste(BuyseSummary$nb$endpoint, " th=",BuyseSummary$nbthreshold ," strata=", BuyseSummary$nb$strata)
    diff <- cbind(diff, diff.sum)
  }
  return(diff)
}

expect_equalPairsBT <- function(BuyseRes1, BuyseRes2){
  count1 <- getCount(BuyseRes1) 
  count2 <- getCount(BuyseRes2)
  expect_equal(count1, count2)
}

expect_equalBT <- function(BuyseRes1, BuyseRes2, slots = NULL, trace = 1){
  if(is.null(slots)){slots <- intersect(names(attributes(BuyseRes1)), names(attributes(BuyseRes2)))}
  
  for(iterSlot in slots){
    res <- try(expect_equal(slot(BuyseRes1,  iterSlot), slot(BuyseRes2,  iterSlot)), silent = TRUE)
    if("try-error" %in% class(res) ){
      cat("Differences in slot: ",iterSlot,"\n")
      if(trace == 2){cat(res[1])}
    }
  }
}