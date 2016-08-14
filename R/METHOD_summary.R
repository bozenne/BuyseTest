#' @docType methods
#' @name BuyseRes-summary
#' @title Summary Method for Class "BuyseRes"
#' @aliases summary BuyseRes-summary
#' @include OBJECT_BuyseTest.R
#' 
#' @description Summarize the results from the \code{\link{BuyseTest}} function.
#' 
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' @param show the type of result to print. Can be \code{"nb"} or \code{"pc"} or \code{NULL} (nothing is printed). Default is \code{"pc"}.
#' @param strata the name of the strata to be displayed \emph{character vector}. Default is \code{"global"} which displays the overall results.
#' @param digit the number of digit to use for printing the results. \emph{integer}. Default is \code{3}. 
#' 
#' @details WARNING : the confidence interval is computed using quantiles of the distribution of the cumulative proportion in favor of the treatment under the null hypothesis. It thus may not be valid if this hypothesis is rejected. In particular, if the cumulative proportion in favor of the treatment is close to 1, the upper limit of the confidence interval may exceed 1. 
#' 
#' @return 
#'   A \code{"List"} composed of two matrices containing the endpoint (and the strata) in rows and the results of the pairwise comparison in columns:
#'     \itemize{
#'       \item\code{[[nb]]} : The favorable, unfavorable, neutral and uninformative pairs are reported in number of pairs classified in each category \cr (\code{"n.favorable","n.unfavorable","n.neutral","n.uninformative"}). 
#'       \item\code{[[pc]]} : The favorable, unfavorable, neutral and uninformative pairs are reported in percentage of pair classified in each category \cr (\code{"pc.favorable","pc.unfavorable","pc.neutral","pc.uninformative"}). 
#'     }
#'   \code{"delta"} indicates for each endpoint (and each strata) thechance of a better outcome and \code{"Delta"} the cumulative chance of a better outcome. \cr
#'   The confidence interval and the p.value are given for the cumulative chance of a better outcome (\code{"CIinf.Delta","CIsup.Delta","p.value"})
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{BuyseRes-class}} for a presentation of the \code{BuyseRes} object.
#' 
#' @examples
#'   n.Treatment_testBin <- 500
#'   n.Control_testBin <- 500
#'   prob.Treatment_testBin <- c(0.5,0.75)
#'   prob.Control_testBin <- c(0.5,0.25)
#'   
#'   set.seed(10)
#'   data_testBin <- data.frame(treatment=c(rep(1,n.Treatment_testBin),rep(0,n.Treatment_testBin)))
#'   data_testBin$endpoint1 <- c(rbinom(n.Treatment_testBin,size=1,prob=prob.Treatment_testBin[1]),
#'                               rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[1]))
#'   data_testBin$endpoint2 <- c(rbinom(n.Control_testBin,size=1,prob=prob.Treatment_testBin[2]),
#'                               rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[2]))
#'   data_testBin$strata <- rbinom(n.Treatment_testBin+n.Control_testBin,size=4,prob=0.5)
#'   
#'   #### no strata
#'   \dontrun{
#'     BuyseTest_object <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
#'                                   treatment="treatment",type,=c("bin","bin"),n.bootstrap=10000)
#'   }
#'   \dontshow{
#'     BuyseTest_object <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
#'                                   treatment="treatment",type=c("bin","bin"),
#'                                   n.bootstrap=10,trace=0)
#'   }
#'   
#'   summary_BuyseTest_object <- summary(BuyseTest_object)
#' 
#' @keywords summary BuyseRes-method

#' @rdname BuyseRes-summary
setGeneric(name = "summary", 
           def = function(object, ...){standardGeneric("summary")}
)


#' @rdname BuyseRes-summary
#' @exportMethod summary
setMethod(f = "summary",
          signature = "BuyseRes",
          definition = function(object, show = "pc", strata = NULL, digit = c(2,3)){
            
            # preparation
            if(!is.null(show) && show %in% c("nb","pc") == FALSE){
              stop("summary[BuyseRes] : wrong specification of \'show\' \n",
                   "valid values : \"nb\" \"pc\" \n",
                   "proposed \'show\' : ",show,"\n")
            }
            
            # mise en forme
            n.endpoint <- length(object@endpoint)
            n.strata <- length(object@strata)
            
            table <- list(nb=data.frame(matrix(NA,nrow=(n.strata+1)*n.endpoint,ncol=15)),
                          pc=data.frame(matrix(NA,nrow=(n.strata+1)*n.endpoint,ncol=15)))
            names(table$nb) <- c("endpoint","threshold","strata","n.total","n.favorable","n.unfavorable","n.neutral","n.uninf","delta","Delta","CIinf.Delta","CIsup.Delta","n.bootstrap","p.value","")
            names(table$pc) <- c("endpoint","threshold","strata","pc.total","pc.favorable","pc.unfavorable","pc.neutral","pc.uninf","delta","Delta","CIinf.Delta","CIsup.Delta","n.bootstrap","p.value","")
            Delta_boot <- apply(apply(object@delta_boot,c(2,3),sum),2,cumsum)
            if(n.endpoint==1){Delta_boot <- rbind(Delta_boot)}
            
            # global
            index.global <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1
            table$nb[index.global,"endpoint"] <- object@endpoint
            table$nb[index.global,"threshold"] <- object@threshold
            table$nb[index.global,"strata"] <- "global"
            table$nb[index.global,"n.favorable"] <- colSums(object@count_favorable)
            table$nb[index.global,"n.unfavorable"] <- colSums(object@count_unfavorable)
            table$nb[index.global,"n.neutral"] <- colSums(object@count_neutral)
            table$nb[index.global,"n.uninf"] <- colSums(object@count_uninf)
            table$nb[index.global,"delta"] <- colSums(object@delta)
            table$nb[index.global,"Delta"] <- cumsum(colSums(object@delta))
            table$nb[index.global,"CIinf.Delta"] <- cumsum(colSums(object@delta))+object@Delta_quantile["2.5%",]
            table$nb[index.global,"CIsup.Delta"] <- cumsum(colSums(object@delta))+object@Delta_quantile["97.5%",]
            table$nb[index.global,"n.bootstrap"] <- apply(Delta_boot,1,function(x){sum(!is.na(x))})
            table$nb[index.global,"p.value"] <- object@p.value
            table$nb[index.global,ncol(table$nb)] <- sapply(object@p.value,function(x){
              if(is.na(x)){NA}else if(x<0.001){"***"}else if(x<0.01){"**"}else if(x<0.05){"*"}else if(x<0.1){"."}else{""}
            })
            table$nb[index.global,"n.total"] <- rowSums(table$nb[index.global,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])
            
            table$pc[index.global,"endpoint"] <- object@endpoint
            table$pc[index.global,"threshold"] <- object@threshold
            table$pc[index.global,"strata"] <- "global"
            table$pc[index.global,"pc.favorable"] <- 100*colSums(object@count_favorable)/table$nb[1,"n.total"]#table$nb[index.global,"n.total"]
            table$pc[index.global,"pc.unfavorable"] <- 100*colSums(object@count_unfavorable)/table$nb[1,"n.total"]#table$nb[index.global,"n.total"]
            table$pc[index.global,"pc.neutral"] <- 100*colSums(object@count_neutral)/table$nb[1,"n.total"]#table$nb[index.global,"n.total"]
            table$pc[index.global,"pc.uninf"] <- 100*colSums(object@count_uninf)/table$nb[1,"n.total"]#table$nb[index.global,"n.total"]
            table$pc[index.global,"delta"] <- colSums(object@delta)
            table$pc[index.global,"Delta"] <- cumsum(colSums(object@delta))
            table$pc[index.global,"CIinf.Delta"] <- cumsum(colSums(object@delta))+object@Delta_quantile["2.5%",]
            table$pc[index.global,"CIsup.Delta"] <- cumsum(colSums(object@delta))+object@Delta_quantile["97.5%",]  
            table$pc[index.global,"p.value"] <- object@p.value
            table$pc[index.global,"n.bootstrap"] <- apply(Delta_boot,1,function(x){sum(!is.na(x))})
            table$pc[index.global,ncol(table$pc)] <- sapply(object@p.value,function(x){
              if(is.na(x)){NA}else if(x<0.001){"***"}else if(x<0.01){"**"}else if(x<0.05){"*"}else if(x<0.1){"."}else{""}
            })
            table$pc[index.global,"pc.total"] <- rowSums(table$pc[index.global,c("pc.favorable","pc.unfavorable","pc.neutral","pc.uninf")]) 
            
            for(iter_strata in 1:n.strata){
              index.strata <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1+iter_strata
              table$nb[index.strata,"endpoint"] <- object@endpoint
              table$nb[index.strata,"threshold"] <- object@threshold
              table$nb[index.strata,"strata"] <- object@strata[iter_strata]
              table$nb[index.strata,"n.favorable"] <- object@count_favorable[iter_strata,]
              table$nb[index.strata,"n.unfavorable"] <- object@count_unfavorable[iter_strata,]
              table$nb[index.strata,"n.neutral"] <- object@count_neutral[iter_strata,]
              table$nb[index.strata,"n.uninf"] <- object@count_uninf[iter_strata,]
              table$nb[index.strata,"delta"] <- object@delta[iter_strata,]
              table$nb[index.strata,"n.total"] <- rowSums(table$nb[index.strata,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])
             
              table$pc[index.strata,"endpoint"] <- object@endpoint
              table$pc[index.strata,"threshold"] <- object@threshold
              table$pc[index.strata,"strata"] <- object@strata[iter_strata]
              table$pc[index.strata,"pc.favorable"] <- 100*object@count_favorable[iter_strata,]/table$nb[1,"n.total"]#table$nb[index.strata,"n.total"]
              table$pc[index.strata,"pc.unfavorable"] <- 100*object@count_unfavorable[iter_strata,]/table$nb[1,"n.total"]#table$nb[index.strata,"n.total"]
              table$pc[index.strata,"pc.neutral"] <- 100*object@count_neutral[iter_strata,]/table$nb[1,"n.total"]#table$nb[index.strata,"n.total"]
              table$pc[index.strata,"pc.uninf"] <- 100*object@count_uninf[iter_strata,]/table$nb[1,"n.total"]#table$nb[index.strata,"n.total"]
              table$pc[index.strata,"delta"] <- object@delta[iter_strata,]      
              table$pc[index.strata,"pc.total"] <- rowSums(table$pc[index.strata,c("pc.favorable","pc.unfavorable","pc.neutral","pc.uninf")])
            }
            
            #### posttreatment
            ## bootstrap
            if(all(table$pc[index.global,"n.bootstrap"]==0)){
              keep.cols <- setdiff(names(table$nb), c("CIinf.Delta","CIsup.Delta","n.bootstrap","p.value",""))
              table$nb <- table$nb[,keep.cols, drop = FALSE]
              keep.cols <- setdiff(names(table$pc), c("CIinf.Delta","CIsup.Delta","n.bootstrap","p.value",""))
              table$pc <- table$pc[,keep.cols, drop = FALSE]
            }
            
            ## strata
            if(!is.null(strata)){
            table$nb <- table$nb[table$nb$strata %in% strata,,drop = FALSE]
            table$pc <- table$pc[table$pc$strata %in% strata,,drop = FALSE]
            }
            
            if(n.strata == 1){
              keep.cols <- which(names(table$nb) %in% setdiff(names(table$nb), "strata"))
              table$nb <- table$nb[,keep.cols,drop = FALSE]
              keep.cols <- which(names(table$pc) %in% setdiff(names(table$pc), "strata"))
              table$pc <- table$pc[,keep.cols,drop = FALSE]
            }
            
            ## rounding
            if(length(digit) == 1){digit <- rep(digit,2)}
            
            if(!is.na(digit[1]) && digit[1]>=0){
              param.signif <- c("pc.total","pc.favorable","pc.unfavorable","pc.neutral","pc.uninf")
              table$pc[,param.signif] <- sapply(table$pc[,param.signif],round,digit=digit[1])
            }
            if("n.bootstrap" %in% names(table$pc) && !is.na(digit[2]) && digit[2]>=0){
              param.signif <- c("delta","Delta","CIinf.Delta","CIsup.Delta")
              table$nb[,param.signif] <- sapply(table$nb[,param.signif],signif,digit=digit[2])
              table$pc[,param.signif] <- sapply(table$pc[,param.signif],signif,digit=digit[2])
            }
            
            # affichage
            if(!is.null(show)){
              table.print <- switch(show,
                                    "nb" = table$nb,
                                    "pc" = table$pc
              )
              emptyCol <- which(names(table.print) %in% c("Delta","CIinf.Delta","CIsup.Delta","n.bootstrap","p.value",""))
              for(iterCol in emptyCol){
                table.print[is.na(table.print[,iterCol]), iterCol] <- ""
              }
              table.print$threshold[is.na(table.print$threshold)] <- ""
              cat("        BuyseTest \n")
              cat("> Groups: ",object@levels.treatment[1],"(control) vs. ",object@levels.treatment[2],"(Treatment) \n")
              if(n.strata>1){cat("> ",n.strata," strata \n", sep = "")}
              cat("> ",n.endpoint," endpoint",if(n.endpoint>1){"s"}, "\n", sep = "")
              print(table.print, row.names = FALSE)         
            }
            
            # export
            return(invisible(table))
            
          }
)
