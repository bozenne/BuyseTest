### S4-BuysePower-model.tables.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2023 (14:29) 
## Version: 
## Last-Updated: jun 28 2023 (13:59) 
##           By: Brice Ozenne
##     Update #: 28
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.tables (documentation)
#' @docType methods
#' @name S4BuysePower-model.tables
#' @title Extract Summary for Class "S4BuysePower"
#' @aliases model.tables,S4BuysePower-method
#' @include S4-BuysePower.R
#' 
#' @description Extract a summary of the results from the \code{\link{powerBuyseTest}} function.
#' 
#' @param x output of \code{\link{powerBuyseTest}}
#' @param statistic [character] statistic relative to which the power should be computed:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016),
#' \code{"mannWhitney"} displays the proportion in favor of the treatment (also called Mann-Whitney parameter), as described in Fay et al. (2018).
#' Default value read from \code{BuyseTest.options()}.
#' @param endpoint [character vector] the endpoints to be displayed: must be the name of the endpoint followed by an underscore and then by the threshold.
#' @param transformation [logical] should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' @param order.Hprojection [integer 1,2] the order of the H-project to be used to compute the variance of the net benefit/win ratio.
#'
#' @seealso 
#' \code{\link{powerBuyseTest}} for performing a simulation study for generalized pairwise comparison. \cr
#'
#' @return data.frame
#' @keywords methods
#' @author Brice Ozenne

## * model.tables (code)
#' @exportMethod model.tables
setMethod(f = "model.tables",
          signature = "S4BuysePower",
          definition = function(x, statistic = NULL, endpoint = NULL, order.Hprojection = NULL, transformation = NULL){

              dt.res <- slot(x, name = "results")
              alpha <- 1-slot(x, name = "conf.level")
              null <- slot(x, name = "null")
              method.inference <- slot(x, name = "method.inference")
              object.endpoint <- slot(x, name = "endpoint")
              object.restriction <- slot(x, name = "restriction")
              object.threshold <- slot(x, name = "threshold")
              object.type <- slot(x, name = "type")

              ## ** normalize and check arguments
              valid.endpoint <- names(object.endpoint)
              valid.statistic <- unique(dt.res$statistic)
              valid.order <- unique(dt.res$order)
              valid.transformation <- unique(dt.res$transformation)
              option <- BuyseTest.options()
              if(is.null(statistic)){
                  statistic <- unique(dt.res$statistic)
              }
              if(is.null(endpoint)){
                  endpoint <- utils::tail(valid.endpoint, 1)
              }else if(identical(endpoint,"all")){
                  endpoint <- valid.endpoint
              }else if(is.numeric(endpoint) && all(endpoint %in% 1:length(valid.endpoint))){
                  endpoint <- valid.endpoint[endpoint]
              }
              if(is.null(order.Hprojection)){
                  order.Hprojection <- max(dt.res$order.Hprojection)
              }
              if(is.null(transformation)){
                  transformation <- any(dt.res$transformation!="none")                  
              }

              statistic <- sapply(gsub("[[:blank:]]", "", tolower(statistic)),
                                  switch,
                                  "netbenefit" = "netBenefit",
                                  "winratio" = "winRatio",
                                  "favorable" = "favorable",
                                  "unfavorable" = "unfavorable",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = valid.statistic,
                             valid.length = 1:2,
                             method = "summary[S4BuysePower]")

              validCharacter(endpoint,
                             name1 = "endpoint",
                             valid.length = NULL,
                             valid.values = valid.endpoint,
                             refuse.duplicates = TRUE,
                             refuse.NULL = TRUE,
                             method = "summary[S4BuysePower]")

              validLogical(transformation,
                           name1 = "transformation",
                           valid.length = 1,
                           method = "summary[S4BuysePower]")

              validInteger(order.Hprojection,
                           name1 = "order.Hprojection",
                           valid.length = 1,
                           min = min(valid.order),
                           max = max(valid.order),
                           method = "summary[S4BuysePower]")

              ## ** subset
              if(transformation){
                  index.subset <- which((dt.res$endpoint %in% endpoint) * (dt.res$order == order.Hprojection) * (dt.res$transformation != "none") == 1)
              }else{
                  index.subset <- which((dt.res$endpoint %in% endpoint) * (dt.res$order == order.Hprojection) * (dt.res$transformation == "none") == 1)
              }
              if(method.inference == "none"){                          
                  dtS.res <- dt.res[index.subset,list(rep.estimate = sum(!is.na(.SD$estimate)),
                                                      mean.estimate = mean(.SD$estimate, na.rm = TRUE)),
                                    by = c("n.T","n.C","endpoint","statistic"),]
                  col.value <- c("mean.estimate","rep.estimate")
              }else{
                  dtS.res <- dt.res[index.subset,list(rep.estimate = sum(!is.na(.SD$estimate)),
                                                      rep.se = sum(!is.na(.SD$se)),
                                                      mean.estimate = mean(.SD$estimate, na.rm = TRUE),
                                                      sd.estimate = stats::sd(.SD$estimate, na.rm = TRUE),
                                                      mean.se = mean(.SD$se, na.rm = TRUE),
                                                      rejection.rate = mean(.SD$p.value<=alpha, na.rm = TRUE)),
                                    by = c("n.T","n.C","endpoint","statistic"),]
                  col.value <- c("mean.estimate","sd.estimate","mean.se","rejection.rate","rep.estimate","rep.se")
              }
              index.endpoint <- match(dtS.res$endpoint, valid.endpoint)
              dtS.res$endpoint <- object.endpoint[index.endpoint]
              dtS.res$threshold <- object.threshold[index.endpoint]
              dtS.res$restriction <- object.restriction[index.endpoint]
              if(any(object.type[index.endpoint]=="bin")){
                  dtS.res$threshold[object.type[index.endpoint]=="bin"] <- NA
              }
              data.table::setkeyv(dtS.res, c("endpoint","n.T"))
              data.table::setcolorder(dtS.res, neworder = c("statistic","endpoint","restriction","threshold","n.T","n.C",col.value))

              ## ** export
              if(method.inference == "u-statistic"){                          
                  if(transformation){
                      attr(dtS.res,"transformation") <- stats::setNames(dt.res[index.subset,.SD$transformation],dt.res[index.subset,.SD$statistic])[!duplicated(dt.res[index.subset,.SD$transformation])]
                  }
                  attr(dtS.res,"order.Hprojection") <- order.Hprojection
              }
              return(dtS.res)
          }
          )

##----------------------------------------------------------------------
### S4-BuysePower-model.tables.R ends here
