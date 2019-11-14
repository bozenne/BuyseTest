### BuyseRes-coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 12 2019 (10:45) 
## Version: 
## Last-Updated: nov 14 2019 (14:41) 
##           By: Brice Ozenne
##     Update #: 38
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - coef
#' @docType methods
#' @name BuyseRes-coef
#' @title Coef Method for Class "BuyseRes"
#' @aliases coef,BuyseRes-method
#' @include BuyseRes-object.R
#' 
#' @description Extract summary statistics from the result of a \code{\link{BuyseTest}} function.
#' 
#' @param object output of \code{\link{BuyseTest}}
#' @param statistic [character] the type of summary statistic. See the detail section.
#' @param stratified [logical] should the summary statistic be strata-specific?
#' Otherwise a summary statistic over all strata is returned.
#' @param cumulative [logical] should the score be cumulated over endpoints?
#' Otherwise display the contribution of each endpoint.
#' @param ... ignored.
#'
#' @details
#' One of the following statistic can be specified:
#' \itemize{
#' \item \code{"netBenefit"}: returns the net benefit.
#' \item \code{"winRatio"}: returns the win ratio.
#' \item \code{"count.favorable"}: returns the number of pairs in favor of the treatment.
#' \item \code{"count.unfavorable"}: returns the number of pairs in favor of the control.
#' \item \code{"count.neutral"}: returns the number of neutral pairs.
#' \item \code{"count.uninf"}: returns the number of uninformative pairs.
#' \item \code{"pc.favorable"}: returns the percentage of pairs in favor of the treatment, i.e. \eqn{P[X \geq Y + \tau]}.
#' \item \code{"pc.unfavorable"}: returns the percentage of pairs in favor of the control, i.e. \eqn{P[Y \geq X + \tau]}.
#' \item \code{"pc.neutral"}: returns the percentage of neutral pairs.
#' \item \code{"pc.uninf"}: returns the percentage of uninformative pairs.
#' }
#' @keywords coef BuyseRes-method
#' @author Brice Ozenne

## * method - coef
#' @rdname BuyseRes-coef
#' @exportMethod coef
setMethod(f = "coef",
          signature = "BuyseRes",
          definition = function(object,
                                statistic = NULL,
                                stratified = FALSE,
                                cumulative = TRUE,
                                ...){

              ## ** normalize arguments
              option <- BuyseTest.options()
              if(is.null(statistic)){
                  statistic <- option$statistic
              }

              type.count <- c("count.favorable","count.unfavorable","count.neutral","count.uninf")
              type.pc <- c("pc.favorable","pc.unfavorable","pc.neutral","pc.uninf")
              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio",type.count,type.pc),
                             valid.length = 1,
                             method = "coef[BuyseRes]")

              ## ** extract information
              if(statistic %in% c("netBenefit","winRatio")){

                  if(stratified){
                      out <- slot(object, paste("delta", statistic , sep = "."))
                  }else{
                      out <- slot(object, paste("Delta", statistic, sep = "."))
                  }
                  
              }else if(statistic %in% type.count){
                  
                  out <- slot(object, statistic)
                  if(cumulative){
                      out <- t(apply(out,1,cumsum))
                  }
                  if(!stratified){
                      out <- colSums(out)
                  }
                  
              }else if(statistic %in% type.pc){

                  out <- slot(object, type.count[match(statistic, type.pc)])/sum(slot(object, "n.pairs"))
                  if(cumulative){
                      out <- t(apply(out,1,cumsum))
                  }
                  if(!stratified){
                      out <- colSums(out)
                  }
                  
              }
              
                 

              return(out)

          })

######################################################################
### BuyseRes-coef.R ends here
