## * Documentation - getCount
#' @name BuyseRes-getCount
#' @title get Method for Class "BuyseRes"
#' @include BuyseRes-object.R
#' @aliases BuyseRes-getCount getCount getCount,BuyseRes-method
#'
#' @description Extract the number of pairs.
#'
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' @param type the type of pairs to be counted. Can be \code{"favorable"}, \code{"unfavorable"}, \code{neutral}, or \code{uninf}. Can also be \code{"all"} to select all of them.
#'
#' @return
#'   A \code{"vector"} containing the number of pairs
#'
#' @keywords getCount BuyseRes-method

## * Method - getCount
#' @rdname BuyseRes-getCount
setMethod(f = "getCount",
          signature = "BuyseRes",
          definition = function(object, type){

            if (missing(type)) {
              type <- c("favorable","unfavorable","neutral","uninf")
            }

            validCharacter(type,
                           valid.length = NULL,
                           valid.values = c("favorable","unfavorable","neutral","uninf"),
                           method = "getCount")

            out <- NULL
            if ("favorable" %in% type) {out <- c(out, favorable = object@count.favorable)}
            if ("unfavorable" %in% type) {out <- c(out, unfavorable = object@count.unfavorable)}
            if ("neutral" %in% type) {out <- c(out, neutral = object@count.neutral)}
            if ("uninf" %in% type) {out <- c(out, uninf = object@count.uninf)}

            return(out)
          }
)
