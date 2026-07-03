### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  3 2026 (14:26) 
## Version: 
## Last-Updated: jul  3 2026 (14:26) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .score.logit
#' @title Score for Logistic Regressions
#' @description Compute the first derivative of the log-likelihood IPCW logistic regressions.
#' @noRd
#'
#' @param object a glm object corresponding to a logistic regression.
#' @param indiv [logical] should the individual contribution be output instead of the total score?
#' 
.score.logit <- function(object, indiv){
    X <- stats::model.matrix(object)
    pi <- stats::predict(object, type = "response")
    Y <- object$y
    W <- object$prior.weights
    out <- .colMultiply_cpp(X, scale = W*(Y - pi))
    colnames(out) <- colnames(X)
    if(indiv){
        return(out)
    }else{
        return(colSums(out))
    }
}


##----------------------------------------------------------------------
### score.R ends here
