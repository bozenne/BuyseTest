### as.data.table.performance.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  9 2021 (10:04) 
## Version: 
## Last-Updated: jun 27 2023 (09:55) 
##           By: Brice Ozenne
##     Update #: 95
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * as.data.table.performance
##' @title Convert Performance Objet to data.table
##' @description Extract the AUC/brier score values or the prediction into a data.table format.
##'
##' @param x object of class \code{"performance"}.
##' @param type [character] either \code{"metric"} to extract AUC/brier score or \code{"prediction"} to extract predictions.
##' @param format [character] should the result be outcome in the long format (\code{"long"}) or in the wide format (\code{"wide"}).
##' Note relevant when using \code{type="metric"}.
##' @param keep.rownames Not used. For compatibility with the generic method.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A data.table object
##' @keywords methods
##' 
##' @export
as.data.table.performance <- function(x, keep.rownames = FALSE, type = "performance", format = NULL, ...){

    ## ** normalize user input
    if(length(type)!=1){
        stop("Argument \'type\' must have length 1.")
    }
    type <- match.arg(type, c("performance",
                              "prediction",paste0("prediction-",names(x$prediction)),
                              "roc",paste0("roc-",names(x$prediction)),
                              "fold"))
    if(!is.null(format)){
        format <- match.arg(format, c("long","wide"))
    }

    ## ** extract data
    if(type=="performance"){
        return(as.data.table(x$performance))
    }else if(type %in% c("prediction","prediction-internal","prediction-external","prediction-cv")){
        x.prediction <- x$prediction
        x.response <- x$response
        out <- NULL
        if(type=="prediction-internal"){
            x.prediction <- x.prediction["internal"]
            x.response <- x.response["internal"]
        }else if(type=="prediction-external"){
            x.prediction <- x.prediction["external"]
            x.response <- x.response["external"]
        }else if(type=="prediction-cv"){
            x.prediction <- x.prediction["cv"]
            x.response <- x.response["cv"]
        }

        for(iType in names(x.prediction)){ ## iType <- names(x.prediction)[3]
            if(iType == "internal"){

                iX.prediction <- data.table(method = "internal", outcome = x.response[[iType]], x.prediction[[iType]], observation = 1:NROW(x.prediction[[iType]]),
                                            repetition = as.numeric(NA), fold = as.numeric(NA))

            }else if(iType == "external"){

                iX.prediction <- data.table(method = "external", outcome = x.response[[iType]], x.prediction[[iType]], observation = 1:NROW(x.prediction[[iType]]),
                                            repetition = as.numeric(NA), fold = as.numeric(NA))

            }else if(iType == "cv"){

                iX <- x.prediction[[iType]]
                iIndex <- attr(iX,"index")
                attr(iX,"index") <- NULL
                iX.prediction <- data.table(method = "cv", outcome = as.numeric(NA),do.call(rbind,lapply(apply(iX,3,list),"[[",1)),do.call(rbind,lapply(apply(iIndex,3,list),"[[",1)))
                iX.prediction$outcome <- x.response[[iType]][iX.prediction$observation]

            }
            out <- rbind(out, iX.prediction)
        }

        if(!is.null(format) && format == "long"){
            out <- data.table::melt(out, id.vars = intersect(names(out),c("method","outcome","observation","repetition","fold")),
                                    variable.name = "model", value.name = "prediction")
        }
        return(out)
    }else if(type %in% c("roc","roc-internal","roc-external","roc-cv")){
        newx <- as.data.table.performance(x, type = gsub("roc","prediction",type), format = "long")
        Umethod <- unique(newx$method)
        Umodel <- unique(x$model)
        out <- NULL
        for(iMethod in Umethod){

            iNewx <- newx[newx$method==iMethod]
            setkeyv(iNewx, c("repetition","model","prediction"))
            ## se: among those who have the outcome P[score>=threshold|Y=1]
            ## sp: among those who do not have the outcome P[score<threshold|Y=0]
            iOut <- iNewx[,list("observation"=c(NA,.SD$observation),
                                "threshold"=c(0,.SD$prediction),
                                "se"=rev(cumsum(c(0,rev(.SD$outcome))==1))/sum(.SD$outcome==1),
                                "sp"=cumsum(c(1,.SD$outcome)==0)/sum(.SD$outcome==0)), ## below threshold classified as 1: sp is the number of 0 divided by the number of negative
                          by = c("repetition","model")]
            out <- rbind(out,iOut)
        }

        if(!is.null(format) && format == "wide"){
            out <- data.table::dcast(out, formula = repetition+observation~model,
                                     value.var = c("threshold","se","sp"),sep=".")
        }else{
            out <- out[order(out[["repetition"]],out$model,1-out$sp,out$se)]
        }
        return(out)
    }else if(type == "fold"){

        if("cv" %in% names(x$prediction) == FALSE){
            message("No fold to extract as no cross-validation was performed. \n")
            return(NULL)
        }

        out <- as.data.table(do.call(rbind,lapply(1:dim(attr(x$prediction$cv,"index"))[3], function(k){attr(x$prediction$cv,"index")[,,k]})))

        return(out)
    }
}

## * as.data.table.performanceResample
##' @export
as.data.table.performanceResample <- function(x, ...){
    return(as.data.table(unclass(x)))

}

##----------------------------------------------------------------------
### as.data.table.performance.R ends here
