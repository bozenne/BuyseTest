## * Documentation - simCompetingRisks
#' @title Simulation of Gompertz competing risks data for the BuyseTest
#' @description Simulate Gompertz competing risks data with proportional (via prespecified sub-distribution hazard ratio) or
#' non-proportional sub-distribution hazards. A treatment variable with two groups (treatment and control) is created.
#' @name simCompetingRisks
#' 
#' @param n.T [integer, >0] number of patients in the treatment arm
#' @param n.C [integer, >0] number of patients in the control arm
#' @param p.1C [integer, >0] proportion of events of interest in the control group. Can be NULL if and only if \code{(b.1T, b.1C, b.2T, b.2C)}
#' are provided.
#' @param sHR [double, >0] pre-specified sub-distribution hazard ratio for event of interest. Can be NULL if and only if 
#' \code{(b.1T, b.1C, b.2T, b.2C)} are provided.
#' @param v.1C,v.1T,v.2C,v.2T [double, <0] shape parameters for Gompertz distribution of time to event of interest in control/treatment (C/T) 
#' group and of time to competing event in control/treatment (C/T) group respectively
#' @param b.1C,b.1T,b.2C,b.2T [double, >0] rate parameters for Gompertz distribution of time to event of interest in control/treatment (C/T) 
#' group and of time to competing event in control/treatment (C/T) group respectively. Can be NULL if and only if \code{(p.1C, sHR)} are 
#' provided.
#' @param cens.distrib [character] censoring distribution. Can be \code{"exponential"} for exponential censoring or \code{"uniform"} for
#' uniform censoring. NULL means no censoring.
#' @param param.cens [>0] parameter for censoring distribution. Should be a double for rate parameter of exponential censoring distribution 
#' or a vector of doubles for lower and upper bounds of uniform censoring distribution. NULL means no censoring
#' @param latent [logical] If \code{TRUE}, also export the latent variables (e.g. true event times, true event types and censoring times). 
#' NULL sets this parameter to \code{FALSE}.
#' 
#' @details 
#' The times to the event of interest and to the competing event in each group follow an improper Gompertz distribution 
#' (see Jeong and Fine, 2006), whose cumulative distribution function is 
#' 
#' F(t; b, v) = 1 - exp(b (1 - exp (v t)) / v) \cr 
#' 
#' and hazard functions is
#' 
#' h(t; b, v) = b exp(v t)\cr 
#' 
#' The shape parameters must be negative to have improper distributions for the times to the two events in each group. Note however that 
#' in each group, the overall cumulative incidence function must be proper (i.e. the maximum values of the cumulative incidence of each 
#' event type sum up to 1 in each group). When only providing the shape parameters, the rate parameters are
#' computed to fulfill this condition. In case you whish to provide the rate parameters too, make sure that the condition is met.
#'
#' @examples
#' 
#' #### Providing p.1C and sHR ####
#' d <- simCompetingRisks(n.T = 100, n.C = 100, p.1C = 0.55, v.1C = -0.30, 
#' v.1T = -0.30, v.2C = -0.30, v.2T = -0.30, sHR = 0.5, b.1T = NULL, 
#' b.1C = NULL, b.2T = NULL, b.2C = NULL)
#' 
#' #### Providing the rate parameters ####
#' d <- simCompetingRisks(n.T = 100, n.C = 100, p.1C = NULL, v.1C = -0.30, 
#' v.1T = -0.30, v.2C = -0.30, v.2T = -0.30, sHR = NULL, b.1T = 0.12, 
#' b.1C = 0.24, b.2T = 0.33, b.2C = 0.18)
#' 
#' #### With exponential censoring ####
#' d <- simCompetingRisks(n.T = 100, n.C = 100, p.1C = 0.55, v.1C = -0.30, 
#' v.1T = -0.30, v.2C = -0.30, v.2T = -0.30, sHR = 0.5, b.1T = NULL, 
#' b.1C = NULL, b.2T = NULL, b.2C = NULL, cens.distrib = "exponential", 
#' param.cens = 0.8, latent = TRUE)
#'
#' ### With uniform censoring ####
#' d <- simCompetingRisks(n.T = 100, n.C = 100, p.1C = 0.55, v.1C = -0.30, 
#' v.1T = -0.30, v.2C = -0.30, v.2T = -0.30, sHR = 0.5, b.1T = NULL, 
#' b.1C = NULL, b.2T = NULL, b.2C = NULL, cens.distrib = "uniform", 
#' param.cens = c(0, 7), latent=TRUE)        
#' 
#' @references Jeong J-H. and Fine J. (2006) \bold{Direct parametric inference for the cumulative incidence function}. \emph{Journal of the Royal Statistical
#' Society} 55: 187-200 \cr
#'
#' @author Eva Cantagallo
#' 
#' @keywords function simulations
#'

## * Function simCompetingRisks
#' @rdname simCompetingRisks
#' @export
#' 
  
simCompetingRisks <- function(n.T, n.C, p.1C = NULL, v.1C, v.1T, v.2C, v.2T, sHR = NULL, 
                              b.1T = NULL, b.1C = NULL, b.2T = NULL, b.2C = NULL,
                              cens.distrib = NULL, param.cens = NULL, latent = NULL) {
  
  # Compute rate parameters if not provided
  if(!is.null(b.1T) & !is.null(b.1C) & !is.null(b.2T) & !is.null(b.2C)) {
    p.1T <- 1 - exp(b.1T / v.1T)
    p.1C <- 1 - exp(b.1C / v.1C)
  } else if(!is.null(p.1C) & !is.null(sHR)) {
    b.1C <- v.1C * log(1 - p.1C)
    b.1T <- b.1C * sHR
    p.1T <- 1 - exp(b.1T / v.1T); p.2T <- 1 - p.1T
    p.2C <- 1 - p.1C
    b.2C <- v.2C * log(1 - p.2C)
    b.2T <- v.2T * log(1 - p.2T)
  } else {
    stop("Missing input argument: please provide either (b.1T, b.1C, b.2T, b.2C) or (p.1C, sHR)")
  }

  rF1T <- function(x) log(1 - v.1T * log(1 - x) / b.1T) / v.1T
  rF1C <- function(x) log(1 - v.1C * log(1 - x) / b.1C) / v.1C
  rF2T <- function(x) log(1 - v.2T * log(1 - x) / b.2T) / v.2T
  rF2C <- function(x) log(1 - v.2C * log(1 - x) / b.2C) / v.2C
  n <- (n.T + n.C)
  u <- stats::runif(n, 0, 1)
  data <- data.frame(treatment = c(rep(1, n.T), rep(0, n.C)), event.time = rep(0, n), event.type = rep(0, n))
  indexT1 <- which(data$treatment == 1 & u < p.1T)
  indexT2 <- which(data$treatment == 1 & u >= p.1T)
  indexC1 <- which(data$treatment == 0 & u < p.1C)
  indexC2 <- which(data$treatment == 0 & u >= p.1C)
  data$event.time[indexT1] <- rF1T(u[indexT1])
  data$event.type[indexT1] <- 1
  data$event.time[indexT2] <- rF2T(u[indexT2] - p.1T)
  data$event.type[indexT2] <- 2
  data$event.time[indexC1] <- rF1C(u[indexC1])
  data$event.type[indexC1] <- 1
  data$event.time[indexC2] <- rF2C(u[indexC2] - p.1C)
  data$event.type[indexC2] <- 2
  
  if(!is.null(cens.distrib)) {
    if(cens.distrib == "exponential") {
      data$censoring.time <- stats::rexp(n, rate = param.cens[1])
    }
    else if (cens.distrib == "uniform") {
      if(is.na(param.cens[2])) {
        stop("Missing parameter for uniform censoring distribution")
      }
      data$censoring.time <- stats::runif(n, min = param.cens[1], max = param.cens[2])
    }
    data$time <- apply(data[, c("event.time", "censoring.time")], 1, min)
    data$status <- ifelse(data$time == data$event.time, data$event.type, 0)
    
    if(!latent | is.null(latent)) {
      data_final <- data[, c('treatment', 'time', 'status')]
    } else {
      data_final <- data
    }
  } else {
    data_final <- data
    colnames(data_final) <- c('treatment', 'time', 'status')
  }
  
  return(data_final)
  
}
