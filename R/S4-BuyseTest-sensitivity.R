### sensitivity.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 31 2021 (14:07) 
## Version: 
## Last-Updated: mar 31 2021 (23:18) 
##           By: Brice Ozenne
##     Update #: 134
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - sensitivity
#' @docType methods
#' @name S4BuyseTest-sensitivity
#' @title  Sensitivity Analysis for the Choice of the Thresholds
#' @aliases sensitivity,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Evaluate the statistic of interest along various thresholds of clinical relevance.
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param threshold [list] a list containing for each endpoint the thresholds to be considered.
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016),
#' \code{"favorable"} displays the proportion in favor of the treatment (also called Mann-Whitney parameter), as described in Fay et al. (2018).
#' \code{"unfavorable"} displays the proportion in favor of the control.
#' Default value read from \code{BuyseTest.options()}.
#' @param null [numeric] right hand side of the null hypothesis (used for the computation of the p-value).
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param transformation [logical]  should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}. Not relevant when using permutations or percentile bootstrap.
#' @examples
#' 
#' \dontrun{
#' ## simulate data
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#'
#' ## with one endpoint
#' ff1 <- treatment ~ TTE(eventtime, status = status, threshold = 0.1)
#' BT1 <- BuyseTest(ff1, data= df.data)
#' se.BT1 <- sensitivity(BT1, threshold = seq(0,2,0.25), band = TRUE)
#' autoplot(se.BT1)
#' 
#' ## with two endpoints
#' ff2 <- update(ff1, .~. + cont(score, threshold = 1))
#' BT2 <- BuyseTest(ff2, data= df.data)
#' se.BT2 <- sensitivity(BT2, threshold = list(eventtime = seq(0,2,0.25), score = 0:2),
#'                       band = TRUE)
#' autoplot(se.BT2)
#' }
#' }


## * Method - sensitivity
#' @rdname S4BuyseTest-sensitivity
#' @exportMethod sensitivity
setMethod(f = "sensitivity",
          signature = "S4BuyseTest",
          definition = function(object, threshold, statistic = NULL, band = FALSE, conf.level = NULL, null = NULL, transformation = NULL, alternative = NULL, adj.p.value = FALSE,
                                method.band = "maxT-integration", seed = NA, ...){

              ## ** normalize user input
              ## band
              if(object@method.inference!="u-statistic"){
                  stop("Cannot compute confidence bands when \'method.inference\' used to obtain the object is not \"u-statistic\". \n")
              }
              
              ## endpoint
              name.endpoint <- object@endpoint
              n.endpoint <- length(name.endpoint)
              option <- BuyseTest.options()
              if(is.null(statistic)){
                  statistic <- option$statistic
              }
              if(is.null(transformation)){
                  transformation <- option$transformation
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
              if(is.matrix(threshold)){
                  if(NCOL(threshold)!=n.endpoint){
                      stop("When a matrix, the argument \'threshold\' should contain ",n.endpoint," columns (and not ",length(threshold),"). \n",
                           "Each column corresponds to a prioritized endpoint. \n")
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
                  
              }else{
                  stop("Argument \'threshold\' should be a list or a matrix \n")
              }

              ## formula
              ls.args <- object@call
              if("formula" %in% names(ls.args)){
                  ls.args$formula <- NULL

                  args.tempo <- initializeFormula(eval(object@call$formula))
                  ls.args[names(args.tempo)] <- args.tempo

              }
              ls.args$trace <- 0

              ## ** run BuyseTest
              grid.threshold <- expand.grid(threshold)
              colnames(grid.threshold) <- name.endpoint
              n.se <- NROW(grid.threshold)
              test.varying <- apply(grid.threshold,2,function(iX){length(unique(iX))!=1})
              if(all(test.varying==FALSE)){
                  stop("Only a single combination of thresholds. No need for a sensitivity analysis.\n")
              }
              gridRed.threshold <- grid.threshold[,which(test.varying),drop=FALSE]
                  
              ls.confint <- vector(mode="list", length = n.se)
              ls.iid <- vector(mode="list", length = n.se)
              
              for(iSe in 1:n.se){
                  iLS.args <- ls.args
                  iLS.args$threshold <- as.double(grid.threshold[iSe,])
                  iBT <- do.call(BuyseTest, args = iLS.args)

                  iConfint <- confint(iBT, statistic = statistic, null = null, conf.level = conf.level, alternative = alternative, transformation = transformation)[n.endpoint,]
                  ls.confint[[iSe]] <- data.frame(c(gridRed.threshold[iSe,,drop=FALSE], iConfint))
                  if(iBT@method.inference=="u-statistic"){
                      ls.iid[[iSe]] <- getIid(iBT, statistic = statistic)
                  }
              }
              
              df.confint <- as.data.frame(do.call(rbind,ls.confint))
              if(object@method.inference=="u-statistic"){
                  attr(df.confint, "iid") <- do.call(cbind,ls.iid)
              }
              
              ## ** compute confidence bands
              if(band){
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
                  if(transformation){
                      type <- switch(statistic,
                                     "netBenefit" = "atanh",
                                     "winRatio" = "log",
                                     "favorable" = "cloglog",## note: note the same transformation as confint
                                     "unfavorable" = "cloglog") ## note: note the same transformation as confint
                  }else{
                      type <- "none"
                  }
                  
                  iBand <- riskRegression::transformCIBP(estimate = rbind(df.confint$estimate),
                                                         se = rbind(df.confint$se),
                                                         iid = A.iid,
                                                         null = null,
                                                         conf.level = conf.level,
                                                         alternative = alternative,
                                                         ci = TRUE, type = type, min.value = min.value, max.value = max.value,
                                                         band = TRUE, method.band = method.band, seed = seed, ...,
                                                         p.value = adj.p.value)

                  attr(df.confint,"quantileBand") <- iBand$quantile
                  df.confint$lower.band <- iBand$lowerBand[1,]
                  df.confint$upper.band <- iBand$upperBand[1,]
                  if(adj.p.value==TRUE){
                      df.confint$adj.p.value <- iBand$adj.p.value[1,]
                  }
              }

              ## ** export
              attr(df.confint,"statistic") <- statistic
              attr(df.confint,"grid") <- grid.threshold
              attr(df.confint,"gridRed") <- gridRed.threshold
              class(df.confint) <- append("sensitivity",class(df.confint))
              return(df.confint)

          })

## * autoplot - sensitivity
autoplot.sensitivity <- function(object, plot = TRUE, col = NULL, ci = TRUE, band = TRUE, label = "Threshold for",
                                 size.line = 1, size.point = 1.75, size.ci = 0.5, alpha = 0.1, ...){

    grid <- attr(object,"gridRed")
    statistic <- switch(attr(object,"statistic"),
                        "netBenefit" = "Net benefit",
                        "winRatio" = "Win ratio",
                        "favorable" = "Proportion of favorable pairs",
                        "unfavorable" = "Proportion of unfavorable pairs")
                        
    if(NCOL(grid)>4){
        stop("No graphical display available when the sensitivity analysis is performed on more than 3 thresholds\n")
    }
    nU.var <- apply(grid,2,function(x){length(unique(x))})
    name.var <- names(sort(nU.var, decreasing = TRUE))
    n.var <- length(name.var)
    if(n.var==1){
        
        if("XXindexXX" %in% names(object)){
            stop("No endpoint should be named \"XXindexXX\" as this name is used internally. \n")
        }
        name.var <- c(name.var,"XXindexXX")
        object <- data.frame(XXindexXX = "1", object)
    }else{
        object[[name.var[2]]] <- factor(object[[name.var[2]]], levels = sort(unique(object[[name.var[2]]])))
    }
    
    ## ** display
    gg <- ggplot2::ggplot(data = object, mapping = ggplot2::aes_string(x = name.var[1], y = "estimate", group = name.var[2]))
    if(band && "lower.band" %in% names(object) && "upper.band" %in% names(object)){
        gg <- gg + ggplot2::geom_ribbon(ggplot2::aes_string(ymin="lower.band", ymax = "upper.band", fill = name.var[2]), alpha = alpha)
    }else{
        band <- FALSE
    }
    if(ci && "lower.ci" %in% names(object) && "upper.ci" %in% names(object)){
        gg <- gg + ggplot2::geom_errorbar(ggplot2::aes_string(ymin="lower.ci", ymax = "upper.ci", color = name.var[2]), size = size.ci)
    }else{
        ci <- FALSE
    }
    gg <- gg + ggplot2::geom_point(ggplot2::aes_string(color = name.var[2]), size = size.point) + ggplot2::geom_line(aes_string(color = name.var[2]), size = size.line)
    gg <- gg + ggplot2::ylab(statistic) + ggplot2::xlab(paste0(label," ",name.var[1]))
    gg <- gg + ggplot2::theme(legend.position = "bottom")
    if(n.var==1){
        if(is.null(col)){
            col <- "black"
        }else if(length(col)!=1){
            stop("Argument \'col\' should have lenght one when the sensitivity analysis is performed on one threshold. \n")
        }
        if(ci){
            gg <- gg + ggplot2::scale_color_manual("CIs", values = col, labels = "")
        }else{
            gg <- gg + ggplot2::scale_color_manual("", values = col, labels = "")
            gg <- gg + ggplot2::guides(color = FALSE)
        }
        if(band){
            gg <- gg + ggplot2::scale_fill_manual("Simulatenous CIs", values = col, labels = "")
        }
        
    }else{
        if(is.null(col)){
            if(ci){
                gg <- gg + ggplot2::labs(color = paste0(label," (CI) ",name.var[2]))
            }else{
                gg <- gg + ggplot2::labs(color = paste0(label," ",name.var[2]))
            }
            if(band){
                gg <- gg + ggplot2::labs(fill = "simulatenous CIs")
            }
        }else{
            if(length(col)!=nU.var[[name.var[2]]]){
                stop("Argument \'col\' should have lenght ",nU.var[[name.var[2]]],", the number of unique thresholds relative to the endpoint \"",name.var[2],"\". \n")
            }
            gg <- gg + ggplot2::scale_color_manual(paste0("threshold of clinical relavance for ",name.var[2]), values = col)
        }
    }

    if(plot){
        print(gg)
    }
    
    return(invisible(gg))
}
##----------------------------------------------------------------------
### sensitivity.R ends here
