### summary.performance.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  6 2022 (14:56) 
## Version: 
## Last-Updated: jun 27 2023 (14:19) 
##           By: Brice Ozenne
##     Update #: 54
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.performance
##' @title Summary Method for Performance Objects
##' @description Summary of the performance of binary classifiers
##'
##' @param object output of performance.
##' @param digits [numeric vector of length 2] number of digits used for the estimates and p-values.
##' @param print [logical] should the performance be printed in the console.
##' @param order.model [character vector] ordering of the models.
##' @param ... not used.
##'
##' @keywords print
##' @method summary performance
##' @export
summary.performance <- function(object, order.model = NULL, digits = c(3,3), print = TRUE, ...){

    ## ** re-order models
    if(!is.null(order.model)){
        if(any(duplicated(order.model))){
            stop("Argument \'order.model\' should not contain duplicated values. \n")
        }
        Umodel <- unique(object$performance$model)
        if(is.numeric(order.model)){
            if(!identical(sort(order.model), 1:length(order.model))){
                stop("Argument \'order.model\' should contain integers from 1 to ",length(order.model)," when numeric. \n")
            }
            order.model <- Umodel[order.model]            
        }else{
            if(any(order.model %in% Umodel == FALSE)){
                stop("Unknown value \"",paste(order.model[order.model %in% Umodel == FALSE], collapse = "\" \""),"\" in argument \'order.model\'. \n")
            }
            if(any(unique(object$performance$model) %in% order.model == FALSE)){
                stop("Missing model \"",paste(Umodel[order.model %in% Umodel == FALSE], collapse = "\" \""),"\" in argument \'order.model\'. \n")
            }
        }
        if(!is.null(object$prediction$internal)){
            object$prediction$internal <- object$prediction$internal[,order.model,drop=FALSE]
        }
        if(!is.null(object$prediction$external)){
            object$prediction$external <- object$prediction$external[,order.model,drop=FALSE]
        }
        if(!is.null(object$prediction$cv)){
            object$prediction$cv <- object$prediction$cv[,order.model,,drop=FALSE]
        }
        if(!is.null(object$resampling)){
            object$resampling[,c("model") := factor(.SD$model,levels = order.model)]
            data.table::setkeyv(object$resampling, c("sample","method","metric","model"))
            object$performance <- .performanceResample_inference(performance = object$performance[order(factor(object$performance$model, levels = order.model)),
                                                                                                  c("method","metric","model","estimate")],
                                                                 resampling = object$resampling,
                                                                 type.resampling = object$args$type.resampling,
                                                                 conf.level = object$args$conf.level)
        }else{
            object$performance <- object$performance[order(factor(object$performance$method, levels = c("internal","external","cv")),
                                                           factor(object$performance$metric, levels = c("auc","brier")),
                                                           factor(object$performance$model,levels=order.model)),,drop=FALSE]
            if(length(order.model)>1){
                object$performance$p.value_comp <- NA

                for(iMethod in unique(object$performance$method)){ ## iMethod <- "internal"
                    for(iMetric in c("auc","brier")){ ## iMetric <- "auc"
                        iIndex <- which(object$performance$method==iMethod & object$performance$metric==iMetric)
                        iBeta <- object$performance[iIndex,"estimate"]
                        iIID <- object[[paste0("iid.",iMetric)]][[iMethod]][,order.model]
                        iStat <- c(NA,diff(iBeta)) / c(NA,sqrt(colSums((iIID[,1:(NCOL(iIID)-1),drop=FALSE]-iIID[,-1,drop=FALSE])^2)))
                        object$performance[iIndex,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                    }
                }
            }
        }
        object$auc <-  lapply(object$auc, function(iL){iL[order.model]})
        object$brier <-  lapply(object$brier, function(iL){iL[order.model]})
    }


    ## ** display results
    df.print <- object$performance

    df.print$p.value <- base::format.pval(df.print$p.value, digits = digits[1], eps = 10^(-digits[2]))
    df.print$p.value[is.na(object$performance$p.value)] <- ""
    df.print$p.value_comp <- base::format.pval(df.print$p.value_comp, digits = digits[1], eps = 10^(-digits[2]))
    df.print$p.value_comp[is.na(object$performance$p.value_comp)] <- ""
    df.print <- df.print[,union(names(which(colSums(!is.na(object$performance))>0)),"estimate")]
    print(df.print, digits = digits[1])

    return(invisible(object$performance))
}

## * summary.performance
##' @method print performance
##' @export
print.performance <- function(x, ...){
    out <- summary(x)
    return(invisible(NULL))
}

##----------------------------------------------------------------------
### summary.performance.R ends here
