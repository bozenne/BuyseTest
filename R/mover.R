### mover.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  7 2025 (12:17) 
## Version: 
## Last-Updated: Jul  7 2025 (13:31) 
##           By: Brice Ozenne
##     Update #: 23
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * mover (documentation)
#' @title Inference using MOVER in match designs (EXPERIMENTAL)
#' @name mover
#'
#' @description Implementation of the method of variance estimates recovery (MOVER), an alternative method for statistic inference in a matched design proposed by Matsouaka et al. (2022), for the Net Treatment Benefit.
#'
#' @param object \R object of class \code{S4BuyseTest}, i.e., output of \code{\link{BuyseTest}}
#' @param endpoint [character] for which endpoint the confidence intervals should be output?
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' @param p.value [logical] should the confidence interval be inverted to obtain a p-value.
#' @param type [character] method used to compute confidence intervals for a proportion: \code{"Wilson"} or \code{"Agresti-Coull"}.
#' @param tol [numeric, >0] numeric tolerance for what is considered neglectable (i.e. 0).
#'
#' @details The p-value calculation is performed by finding the confidence level such that the upper or lower bound of the corresponding confidence interval is 0.
#' This was not part of the methodology presented in the original paper.
#'
#' @references R A Matsouaka, A Coles (2022). \bold{Robust statistical inference for the matched net benefit and the matched win ratio using prioritized composite endpoints}. \emph{Stat Methods Med Res} 31(8):1423-1438. doi: 10.1177/09622802221090761.

## * mover (code)
mover <- function(object, endpoint = NULL, conf.level = 0.95, p.value = TRUE, type = "Wilson", tol = 1e-6){
    
    ## ** normalize user input
    q.low <- stats::qnorm((1-conf.level)/2)
    q.high <- stats::qnorm(1-(1-conf.level)/2)
    type <- match.arg(type, c("Wilson","Agresti-Coull"))
    
    valid.endpoint <- object@endpoint
    if(is.null(endpoint)){
        endpoint <- utils::tail(valid.endpoint,1)
    }else{
        if(length(endpoint)>1){
            stop("Argument \'endpoint\' should have length 1. \n")
        }
        if(endpoint %in% valid.endpoint == FALSE){
            stop("Invalid value for argument \'endpoint\'. \n",
                 "Should be one of: \"",paste(valid.endpoint, collapse = "\", \""),"\". \n")
        }
    }

    ## ** normalize user input
    N <- nobs(object)["pairs"]
    pi_w <- unname(coef(object, endpoint = endpoint, statistic = "favorable"))
    pi_l <- unname(coef(object, endpoint = endpoint, statistic = "unfavorable"))

    if(conf.level == 1){
        DeltaL <- -1
        DeltaU <- 1
    }else{
        Ntilde <- N + q.high^2
        pitilde_w <- (N*pi_w+0.5*q.high^2)/Ntilde
        pitilde_l <- (N*pi_l+0.5*q.high^2)/Ntilde

        if(any(c(pi_w,pi_l)<tol) || any(c(pi_w,pi_l)>1-tol)){
            cor_wl <- 0
        }else{
            cor_wl <- - pi_w * pi_l / sqrt(pi_w * (1 - pi_w) * pi_l * (1 - pi_l))
        }

        ## Wilson
        if(type == "Wilson"){
            L_w <- pitilde_w + 0.5/Ntilde * q.low * sqrt(q.low^2 + 4*N*pi_w*(1-pi_w))
            U_w <- pitilde_w + 0.5/Ntilde * q.high * sqrt(q.high^2 + 4*N*pi_w*(1-pi_w))
            L_l <- pitilde_l + 0.5/Ntilde * q.low * sqrt(q.low^2 + 4*N*pi_l*(1-pi_l))
            U_l <- pitilde_l + 0.5/Ntilde * q.high * sqrt(q.high^2 + 4*N*pi_l*(1-pi_l))
        }else if(type == "Agresti-Coull"){
            L_w <- pitilde_w + q.low * sqrt(pitilde_w*(1-pitilde_w)/Ntilde)
            U_w <- pitilde_w + q.high * sqrt(pitilde_w*(1-pitilde_w)/Ntilde)
            L_l <- pitilde_l + q.low * sqrt(pitilde_l*(1-pitilde_l)/Ntilde)
            U_l <- pitilde_l + q.high * sqrt(pitilde_l*(1-pitilde_l)/Ntilde)
        }

        DeltaL <- (pi_w-pi_l) - sqrt((pi_w - L_w)^2 + (U_l - pi_l)^2 - 2 * cor_wl * (pi_w - L_w) * (U_l - pi_l))
        DeltaU <- (pi_w-pi_l) + sqrt((U_w - pi_w)^2 + (pi_l - L_l)^2 - 2 * cor_wl * (U_w - pi_w) * (pi_l - L_l))
    }

    ## p-value
    if(p.value){
        if(abs(DeltaL)<tol | abs(DeltaU)<tol){ ## CI touch 0
            p <- 1-conf.level
        }else if(DeltaL>0 && DeltaU>0){ ## CI do not overlap 0: p>level
            res.search <- stats::uniroot(f = function(mylevel){mover(object, conf.level = mylevel, p.value = FALSE)["lower"]},
                                  lower = conf.level, upper = 1)
        }else if(DeltaL<0 && DeltaU<0){ ## CI do not overlap 0: p>level
            res.search <- stats::uniroot(f = function(mylevel){mover(object, conf.level = mylevel, p.value = FALSE)["upper"]},
                                  lower = conf.level, upper = 1)
        }else if(DeltaL<0 && DeltaU>0){ ## CI overlaps 0: p<level
            res.search <- stats::uniroot(f = function(mylevel){
                tempo <- mover(object, conf.level = mylevel, p.value = FALSE)[c("lower","upper")]
                return(tempo[which.min(abs(tempo))])
            }, lower = 0, upper = conf.level)
        }
        p.value <- 1-res.search$root            
    }else{
        p.value <- NA
    }

    ## export
    out <- c(estimate = pi_w-pi_l, lower = unname(DeltaL), upper = unname(DeltaU), pvalue = p.value)
    return(out)
}


##----------------------------------------------------------------------
### mover.R ends here
