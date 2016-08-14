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
  if(is.null(slots)){slots <- setdiff(intersect(names(attributes(BuyseRes1)), names(attributes(BuyseRes2))), "class")}
  
  test.strata1 <- !("try-error" %in% class(try(BuyseRes1@strata, silent = TRUE)))
  test.strata2 <- !("try-error" %in% class(try(BuyseRes2@strata, silent = TRUE)))

  if (test.strata1 && length(BuyseRes1@strata) > 1) {
    if(test.strata2 == FALSE){
      
      if(!identical(BuyseRes1@delta, BuyseRes2@delta)){
        if(trace>0){cat("* reorder according to strata in BuyseRes 1 and 2 \n")}
        index.match <- match(BuyseRes1@delta[,1], BuyseRes2@delta[,1]) # could be wrong to consider only the first outcome
        BuyseRes2@delta <- BuyseRes2@delta[index.match,,drop=FALSE]
        BuyseRes2@count_favorable <- BuyseRes2@count_favorable[index.match,,drop=FALSE]
        BuyseRes2@count_unfavorable <- BuyseRes2@count_unfavorable[index.match,,drop=FALSE]
        BuyseRes2@count_neutral <- BuyseRes2@count_neutral[index.match,,drop=FALSE]
        BuyseRes2@count_uninf <- BuyseRes2@count_uninf[index.match,,drop=FALSE]
        BuyseRes2@delta_boot <- BuyseRes2@delta_boot[index.match,,,drop = FALSE]
        
        BuyseRes2@index_neutralT <- sort(BuyseRes2@index_neutralT)
        BuyseRes2@index_neutralC <- sort(BuyseRes2@index_neutralC)
        BuyseRes2@index_uninfT <- sort(BuyseRes2@index_uninfT)
        BuyseRes2@index_uninfC <- sort(BuyseRes2@index_uninfC)
        
        BuyseRes1@index_neutralT <- sort(BuyseRes1@index_neutralT)
        BuyseRes1@index_neutralC <- sort(BuyseRes1@index_neutralC)
        BuyseRes1@index_uninfT <- sort(BuyseRes1@index_uninfT)
        BuyseRes1@index_uninfC <- sort(BuyseRes1@index_uninfC)
     }
    }
    
  }
  
  test.error <- FALSE
  for(iterSlot in slots){
    res <- try(expect_equal(slot(BuyseRes1,  iterSlot), slot(BuyseRes2,  iterSlot)), silent = TRUE)
    if("try-error" %in% class(res) ){
      test.error <- TRUE
      cat("Differences in slot: ",iterSlot,"\n")
      if(trace == 2){cat(res[1])}
    }
  }
  if(test.error){stop("difference between the slots of the two BuyseRes objects \n")}
}


Vexpect_less_than <- function(x,y,...){
  sapply(x, function(X){expect_less_than(X,y,...)})
  return(invisible(TRUE))
}
Vexpect_more_than <- function(x,y,...){
  sapply(x, function(X){expect_more_than(X,y,...)})
  return(invisible(TRUE))
}
Vexpect_equal <- function(x,y,...){
  sapply(x, function(X){expect_equal(X,y,...)})
  return(invisible(TRUE))
}
Vexpect_NA <- function(x,...){
  sapply(x, function(X){expect_true(is.na(X),...)})
  return(invisible(TRUE))
}