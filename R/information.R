### information.R --- 
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

## * .information.wglm
#' @title Information for Logistic Regressions
#' @description Compute the information (i.e. opposit of the expectation of the second derivative of the log-likelihood) for logistic regressions.
#' @noRd
#'
#' @param object a glm object corresponding to a logistic regression.
#' @param indiv [logical] should the individual contribution be output instead of the total information?
#' @param center [logical] should the individual contribution be centered around the average?
#' 
.information.logit <- function(object, indiv, center){
    X <- stats::model.matrix(object)
    n.obs <- NROW(X)
    n.param <- NCOL(X)
    pi <- stats::predict(object, type = "response")
    W <- object$prior.weights
    tXWpi <- t(.colMultiply_cpp(X, scale = W*pi*(1-pi)))

    if(indiv){
        if(center){
            Info <- tXWpi %*% X
            out <- array(unlist(lapply(1:n.obs, function(i){tXWpi[,i,drop=FALSE] %*% X[i,,drop=FALSE]-Info/n.obs})),
                         dim = c(n.param,n.param,n.obs), dimnames = list(colnames(X),colnames(X),NULL))
        }else{
            out <- array(unlist(lapply(1:n.obs, function(i){tXWpi[,i,drop=FALSE] %*% X[i,,drop=FALSE]})),
                         dim = c(n.param,n.param,n.obs), dimnames = list(colnames(X),colnames(X),NULL))
        }
        return(out)
        ## apply(out, MARGIN = 1:2, FUN = sum) - tXWpi %*% X
        
    }else{
        out <- tXWpi %*% X
        rownames(out) <- colnames(out)
        return(out)
    }
}



##----------------------------------------------------------------------
### information.R ends here
