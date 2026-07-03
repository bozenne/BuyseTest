### sensitivity.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 31 2021 (14:07) 
## Version: 
## Last-Updated: jul 23 2025 (16:54) 
##           By: Brice Ozenne
##     Update #: 361
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * sensitivity (documentation)
#' @docType methods
#' @name sensitivity
#' @title Sensitivity Analysis for the Choice of the Thresholds
#' @aliases sensitivity,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Evaluate a summary statistic (net benefit, win ratio, ...) using GPC along various thresholds of clinical relevance.
#' 
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param threshold [list] a list containing for each endpoint the thresholds to be considered.
#' @param statistic [character] the statistic summarizing the pairwise comparison: \code{"netBenefit"}, \code{"winRatio"}, \code{"favorable"}, \code{"unfavorable"}.
#' See the documentation of the \code{coef} method for further details.
#' Default value read from \code{BuyseTest.options()}. 
#' @param null [numeric] right hand side of the null hypothesis (used for the computation of the p-value).
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param band [logical] should simulateneous confidence intervals be computed?
#' @param adj.p.value [logical] should p-value adjusted for multiple comparisons be computed?
#' @param transformation [logical]  should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}. Not relevant when using permutations or percentile bootstrap.
#' @param cpus [integer, >0] the number of CPU to use. Default value is 1.
#' @param trace [logical] Should the execution of the function be traced?
#' @param ... argument passsed to the function \code{transformCIBP} of the riskRegression package.
#'
#' @details Simulateneous confidence intervals and adjusted p-values are computed using a single-step max-test approach via the function \code{transformCIBP} of the riskRegression package.
#'
#' @return An S3 object of class \code{S3sensitivity}.
#' @keywords htest
#' 

## * sensitivity (example)
#' @examples
#' 
#' \dontrun{
#' require(ggplot2)
#' 
#' ## simulate data
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#'
#' ## with one endpoint
#' ff1 <- treatment ~ TTE(eventtime, status = status, threshold = 0.1)
#' BT1 <- BuyseTest(ff1, data= df.data)
#' se.BT1 <- sensitivity(BT1, threshold = seq(0,2,0.25), band = TRUE)
#' plot(se.BT1)
#'
#' ## with two endpoints
#' ff2 <- update(ff1, .~. + cont(score, threshold = 1))
#' BT2 <- BuyseTest(ff2, data= df.data)
#' se.BT2 <- sensitivity(BT2, threshold = list(eventtime = seq(0,2,0.25), score = 0:2),
#'                       band = TRUE)
#' plot(se.BT2)
#' plot(se.BT2, col = NA)
#' }


## * sensitivity (method)
#' @rdname sensitivity
#' @exportMethod sensitivity
setMethod(f = "sensitivity",
          signature = "S4BuyseTest",
          definition = function(object, threshold,
                                statistic = NULL, band = FALSE, conf.level = NULL, null = NULL, transformation = NULL, alternative = NULL, adj.p.value = FALSE,
                                trace = TRUE, cpus = 1, ...){

              ## ** normalize user input
              ## band
              if(object@method.inference!="u statistic"){
                  stop("Cannot compute confidence bands when \'method.inference\' used to obtain the object is not \"u-statistic\". \n")
              }
              
              ## endpoint
              name.endpoint <- object@endpoint
              n.endpoint <- length(name.endpoint)
              option <- BuyseTest.options()
              if(is.null(statistic)){
                  statistic <- option$statistic
              }
              if(is.null(conf.level)){
                  conf.level <- option$conf.level
              }
              if(is.null(alternative)){
                  alternative <- option$alternative
              }
              if(is.null(null)){
                  null <- switch(statistic,
                                 "netBenefit" = 0,
                                 "winRatio" = 1,
                                 "favorable" = 1/2,
                                 "unfavorable" = 1/2)
              }else{
                  validNumeric(null, valid.length = 1,
                               min = if("statistic"=="netBenefit"){-1}else{0},
                               max = if("statistic"=="winRatio"){Inf}else{1})
              }

              ## threshold
              if(is.matrix(threshold) || is.data.frame(threshold)){
                  if(is.matrix(threshold)){
                      threshold <- as.data.frame(threshold)
                  }
                  if(NCOL(threshold)!=n.endpoint){
                      stop("When a matrix, the argument \'threshold\' should contain ",n.endpoint," columns (and not ",length(threshold),"). \n",
                           "Each column corresponds to a prioritized endpoint. \n")
                  }
                  grid.threshold <- threshold
                  if(is.null(names(grid.threshold))){
                      names(grid.threshold) <- names(name.endpoint)
                  }else{
                      if(any(duplicated(names(grid.threshold)))){
                          stop("Duplicated column names in argument \"threshold\": \"",paste0(names(grid.threshold)[duplicated(names(grid.threshold))], collapse= "\" \""),"\".\n")
                      }
                      if(any(names(grid.threshold) %in% names(name.endpoint) == FALSE)){
                          stop("Incorrect column names in argument \"threshold\": \"",paste0(names(grid.threshold)[names(grid.threshold) %in% names(name.endpoint) == FALSE], collapse= "\" \""),"\".\n",
                               "Valid names: \"",paste0(setdiff(names(name.endpoint), names(grid.threshold)), collapse= "\" \""),"\".\n")
                      }
                      grid.threshold <- grid.threshold[,names(name.endpoint)]
                  }
                  
              }else if(is.list(threshold) || is.vector(threshold)){

                  if(any(duplicated(name.endpoint))){
                      stop("Argument \'threshold\' must be a matrix when some endpoints are repeteadly used (i.e at different priorities) with different thresholds. \n")
                  }

                  if(!is.list(threshold)){
                      threshold <- list(threshold)
                  }

                  if(is.null(names(threshold))){
                      if(length(threshold)!=n.endpoint){
                          stop("When a list, the argument \'threshold\' should have length ",n.endpoint," (and not ",length(threshold),"). \n",
                               "Each element of the list corresponds to a prioritized endpoint. \n")
                      }
                  }else{
                      if(any(duplicated(names(threshold)))){
                          stop("Argument \'threshold\' must not contain duplicated names. \n",
                               "Duplicated names: \"",paste0(names(threshold)[duplicated(names(threshold))], collapse = "\" \""),"\". \n")
                      }
                      if(any(names(threshold) %in% name.endpoint == FALSE)){
                          stop("Some names used in the argument \'threshold\' does not match the existing endpoints. \n",
                               "Incorrect names: \"",paste0(names(threshold)[names(threshold) %in% name.endpoint == FALSE], collapse = "\" \""),"\". \n",
                               "Possible names: \"",paste0(setdiff(name.endpoint,names(threshold)), collapse = "\" \""),"\". \n")
                      }
                      threshold.save <- threshold
                      threshold <- setNames(vector(mode = "list", length = n.endpoint), name.endpoint)
                      threshold[names(threshold.save)] <- threshold.save
                  }

                  if(any(sapply(threshold,length)==0)){
                      threshold[sapply(threshold,length)==0] <- object@threshold[sapply(threshold,length)==0]
                  }
              
                  grid.threshold <- expand.grid(threshold)
                  colnames(grid.threshold) <- name.endpoint
              
              }else{
                  stop("Argument \'threshold\' should be a list or a matrix \n")
              }

              ## formula
              ls.args <- object@call
              if("formula" %in% names(ls.args)){
                  ls.args$formula <- NULL

                  args.tempo <- initializeFormula(object@call$formula, hierarchical = object@hierarchical)
                  ls.args[setdiff(names(args.tempo),"match")] <- args.tempo[setdiff(names(args.tempo),"match")]
                  if(!is.null(ls.args$strata)){
                      attr(ls.args$strata,"match") <- args.tempo$match
                  }
                  

              }
              ls.args$trace <- 0

              if (cpus == "all") { 
                  cpus <- parallel::detectCores() # this function detect the number of CPU cores 
              }

              if(band && any(object@weightObs!=1)){
                  stop("Confidence bands cannot not currently be derived with weighted observations. \n") 
              }

              ## ** run BuyseTest
              n.se <- NROW(grid.threshold)
              test.varying <- apply(grid.threshold,2,function(iX){length(unique(iX))!=1})
              if(all(test.varying==FALSE)){
                  stop("Only a single combination of thresholds. No need for a sensitivity analysis.\n")
              }
              gridRed.threshold <- grid.threshold[,which(test.varying),drop=FALSE]
                  

              if(trace>0){cat("Run ",n.se," GPC analyses: \n", sep = "")}
              
              if (cpus == 1) {
                  ls.confint <- vector(mode="list", length = n.se)
                  ls.iid <- vector(mode="list", length = n.se)

                  if(trace>0){pb <- utils::txtProgressBar(max = n.se, style = 3)}

                  for(iSe in 1:n.se){
                      if(trace>0){utils::setTxtProgressBar(pb, iSe)}
                      iLS.args <- ls.args
                      iLS.args$threshold <- as.double(grid.threshold[iSe,])
                      iBT <- do.call(BuyseTest, args = iLS.args)

                      iConfint <- confint(iBT, statistic = statistic, null = null, conf.level = conf.level, alternative = alternative, transformation = transformation)[n.endpoint,]
                      ls.confint[[iSe]] <- data.frame(c(gridRed.threshold[iSe,,drop=FALSE], iConfint))
                      if(iBT@method.inference=="u statistic"){
                          ls.iid[[iSe]] <- getIid(iBT, statistic = statistic,simplify=FALSE)$global[,n.endpoint]
                      }
                  }

                  if(trace>0){close(pb)}

              }else{
                  if(trace>0){
                      cl <- suppressMessages(parallel::makeCluster(cpus, outfile = ""))
                      pb <- utils::txtProgressBar(max = n.se, style = 3)          
                  }else{
                      cl <- parallel::makeCluster(cpus)
                  }
                  test.lazyeval <- sapply(ls.args,function(x){inherits(x,"name")})
                  if(any(test.lazyeval)){
                      toExport <- unlist(lapply(ls.args[test.lazyeval],deparse))
                  }else{
                      toExport <- NULL
                  }

                  i <- NULL ## [:forCRANcheck:] foreach
                  ls.sensitivity <- foreach::`%dopar%`(
                                                 foreach::foreach(i=1:n.se, .export = toExport), {                                           
                                                     if(trace>0){utils::setTxtProgressBar(pb, i)}
                                                     iLS.args <- ls.args
                                                     iLS.args$threshold <- as.double(grid.threshold[i,])
                                                     iBT <- do.call(BuyseTest, args = iLS.args)

                                                     iConfint <- confint(iBT, statistic = statistic, null = null, conf.level = conf.level, alternative = alternative, transformation = transformation)[n.endpoint,]
                                                     iOut <- list(confint = data.frame(c(gridRed.threshold[i,,drop=FALSE], iConfint)))
                                                     if(iBT@method.inference=="u statistic"){
                                                         iOut[["iid"]] <- getIid(iBT, statistic = statistic)[,n.endpoint]
                                                     }
                                                     return(iOut)
                                                 })

                  parallel::stopCluster(cl)
                  if(trace>0){close(pb)}

                  ls.confint <- lapply(ls.sensitivity,"[[","confint")
                  if(object@method.inference=="u statistic"){
                      ls.iid <- lapply(ls.sensitivity,"[[","iid")
                  }
              }
              
              df.confint <- as.data.frame(do.call(rbind,ls.confint))
              if(object@method.inference=="u statistic"){
                  attr(df.confint, "iid") <- do.call(cbind,ls.iid)
              }
              
              ## ** compute confidence bands
              if(band || adj.p.value){
                  requireNamespace("riskRegression")
                  A.iid <- array(NA, dim = c(NROW(attr(df.confint, "iid")), NCOL(attr(df.confint, "iid")),1))
                  A.iid[,,1] <- attr(df.confint, "iid")

                  min.value <- switch(statistic,
                                      "netBenefit" = -1,
                                      "winRatio" = 0,
                                      "favorable" = 0,
                                      "unfavorable" = 0)
                  max.value <- switch(statistic,
                                      "netBenefit" = 1,
                                      "winRatio" = Inf,
                                      "favorable" = 1,
                                      "unfavorable" = 1)

                  ## temporary fix: the next few lines should be remove when riskRegression will be updated
                  if(is.null(transformation) || identical(transformation,TRUE)){
                      if(statistic %in% c("none","netBenefit","winRatio") || !inherits(try(riskRegression::transformCIBP(estimate = 1, se = 1, type = "atanh2", seed = NA, band = FALSE, alternative = "two.sided"),silent=TRUE),"try-error")){
                          type <- switch(statistic,
                                         "netBenefit" = "atanh",
                                         "winRatio" = "log",
                                         "favorable" = "cloglog",## note: not the same transformation as confint
                                         "unfavorable" = "cloglog", ## note: not the same transformation as confint
                                         "none" = "none") 
                      }else{
                          type <- switch(statistic,
                                         "netBenefit" = "atanh",
                                         "winRatio" = "log",
                                         "favorable" = "atanh2",
                                         "unfavorable" = "atanh2",
                                         "none" = "none") 
                      }
                  }else if(identical(transformation,FALSE)){
                      type <- "none"
                  }else{
                      type <- transformation                      
                  }
                  dots <- list(...)
                  if("seed" %in% names(dots) == FALSE){
                      dots$seed <- NA
                  }
                  if("method.band" %in% names(dots) == FALSE){
                      dots$method.band <- "maxT-integration"
                  }
                  if("n.sim" %in% names(dots) == FALSE){
                      dots$n.sim <- 10^4
                  }
                  if(any(df.confint$se>0)){
                      iBand <- do.call(riskRegression::transformCIBP,
                                       args = c(list(estimate = rbind(df.confint$estimate[df.confint$se>0]),
                                                     se = rbind(df.confint$se[df.confint$se>0]),
                                                     iid = A.iid[,df.confint$se>0,,drop=FALSE],
                                                     null = null,
                                                     conf.level = conf.level,
                                                     alternative = alternative,
                                                     ci = TRUE, type = type, min.value = min.value, max.value = max.value,
                                                     band = TRUE, p.value = adj.p.value),
                                                dots))
                      if(band){
                          attr(df.confint,"quantileBand") <- iBand$quantile
                          df.confint$lower.band <- rep(0,length(df.confint$se))
                          df.confint$lower.band[df.confint$se>0] <- iBand$lowerBand[1,]
                          df.confint$upper.band <- rep(0,length(df.confint$se))
                          df.confint$upper.band[df.confint$se>0] <- iBand$upperBand[1,]
                      }
                      if(adj.p.value==TRUE){
                          df.confint$adj.p.value <- rep(1,length(df.confint$se))
                          df.confint$adj.p.value[df.confint$se>0] <- iBand$adj.p.value[1,]
                      }
                  }else{
                      if(band && adj.p.value){
                          txt.warning <- "adjusted p-values/confidence bands"
                      }else if(band){
                          txt.warning <- "confidence bands"
                      }else if(band){
                          txt.warning <- "adjusted p-values"
                      }
                      warning("Could not evaluate txt.warning as all standard were estimated to be 0. \n")
                  }

              }

              ## ** export
              attr(df.confint,"statistic") <- statistic
              attr(df.confint,"grid") <- grid.threshold
              attr(df.confint,"gridRed") <- gridRed.threshold
              class(df.confint) <- append("S3sensitivity",class(df.confint))
              return(df.confint)

          })


##----------------------------------------------------------------------
### sensitivity.R ends here
