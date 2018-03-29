## * verboseContext
#' @description Modified context that also display its desc argument
verboseContext <- function(desc){
  context(desc)
  cat(desc,"\n")
  return(invisible(TRUE))
}

## * expect_equalBT
#' @description Test the equality of certain slots between two BuyseTest objects. Usefull when the definition of the object has changed
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
    slot1 <- slot(BuyseRes1,  iterSlot)
    slot2 <- slot(BuyseRes2,  iterSlot)
    
    ## convert from old version to new one
    if(is.list(slot1) && !is.list(slot2)){
      slot1 <- slot1$netChance 
    }else if(is.list(slot1) && !is.list(slot2)){
      slot2 <- slot2$netChance 
    }
    if(iterSlot == "threshold"){
      slot1[is.na(slot1)] <- 0.5
      slot2[is.na(slot2)] <- 0.5
    }
    
    ## test
    res <- try(expect_equal(slot1, slot2, label = iterSlot), silent = TRUE)
    if("try-error" %in% class(res) ){
      test.error <- TRUE
      cat("Differences in slot: ",iterSlot,"\n")
      if(trace == 2){cat(res[1])}
    }
  }
  if(test.error){stop("difference between the slots of the two BuyseRes objects \n")}
}

## * expect_equalPairsBT
#' @description Test the equality of the number of pairs found by two BuyseTest objects
expect_equalPairsBT <- function(BuyseRes1, BuyseRes2){
  count1 <- getCount(BuyseRes1) 
  count2 <- getCount(BuyseRes2)
  expect_equal(count1, count2)
}


## * validPairs
#' @description Test whether the number of pairs found by the summary function is consistent.
#' i.e. the sum over in favor, in defavor, neutral and non informative matches the total number of possible pairs 
#'      the same number of pairs if founded for all strata (expected if the input data is the same for all strata)
validPairs <- function(BuyseRes, type = c("strata","sum")){
  
  BuyseSummary <- summary(BuyseRes, show = FALSE, percentage = FALSE, digit = rep(NA,2))
  enpoint_threshold <- paste(BuyseSummary$endpoint,BuyseSummary$threshold, sep = "_")
  endpoints <- unique(enpoint_threshold)
  D <- length(endpoints)
  index_strata <- which(BuyseSummary$strata!="global")
  levels.strata <- unique(BuyseSummary$strata[index_strata])
  n.strata <- length(levels.strata) 
  
  diff <- NULL
  if("strata" %in% type){
    diff.strata <- matrix(NA,nrow=(1+n.strata)*D,ncol=5)
    colnames(diff.strata) <- c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")
    rownames(diff.strata) <- paste(BuyseSummary$endpoint, " th=",BuyseSummary$threshold ," strata=", BuyseSummary$strata)
    for(iter_endpoint in 1:D){
      index_endpoint <- which(enpoint_threshold==endpoints[iter_endpoint])  
      res_strata <- as.matrix(BuyseSummary[intersect(index_strata,index_endpoint),c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")])
      diff.strata[intersect(index_strata,index_endpoint),] <- t(apply(res_strata,1,function(x){x-res_strata[1,,drop=TRUE]}))
    }
    diff <- cbind(diff, diff.strata)
  }
  if("sum" %in% type){
    diff.sum <- cbind(sum = BuyseSummary$n.total - (BuyseSummary$n.favorable+BuyseSummary$n.unfavorable+BuyseSummary$n.neutral+BuyseSummary$n.uninf))
    rownames(diff.sum) <- paste(BuyseSummary$endpoint, " th=",BuyseSummary$threshold ," strata=", BuyseSummary$strata)
    diff <- cbind(diff, diff.sum)
  }
  return(diff)
}

## * validDelta
#' @description Test whether the number of pairs found by the summary function is consistent with the displayed delta.
validDelta <- function(BuyseRes){
  
  BuyseRes@delta$netChance - (BuyseRes@count_favorable-BuyseRes@count_unfavorable)/BuyseRes@n_pairs
  BuyseRes@Delta$netChance - cumsum(colSums(BuyseRes@delta$netChance))
  
  BuyseRes@delta$winRatio - BuyseRes@count_favorable / (BuyseRes@count_favorable + BuyseRes@count_unfavorable)
  BuyseRes@Delta$winRatio - cumsum(colSums(BuyseRes@count_favorable)) / cumsum(colSums(BuyseRes@count_favorable + BuyseRes@count_unfavorable))
  
  }
  

## * Vexpect_less_than
#' @description Vectorial version of testthat functions
Vexpect_less_than <- function(x,y,...){
  sapply(x, function(X){expect_less_than(X,y,...)})
  return(invisible(TRUE))
}
## * Vexpect_more_than
Vexpect_more_than <- function(x,y,...){
  sapply(x, function(X){expect_more_than(X,y,...)})
  return(invisible(TRUE))
}
## * Vexpect_equal
Vexpect_equal <- function(x,y,...){
  sapply(x, function(X){expect_equal(X,y,...)})
  return(invisible(TRUE))
}
## * Vexpect_NA
Vexpect_NA <- function(x,...){
  sapply(x, function(X){expect_true(is.na(X),...)})
  return(invisible(TRUE))
}
