## * Documentation - simCompetingRisks
#' @name Simulation function
#' @rdname simulation_CR 
#' @title Simulation of Gompertz competing risks data for the BuyseTest
#' @aliases simCompetingRisks
#' 
#' @description Simulate Gompertz competing risks data with proportional (via prespecified sub-distribution hazard ratio) or
#' non-proportional sub-distribution hazards. A treatment variable with two groups (treatment and control) is created. 
#' @param n.T [integer, >0] number of patients in the treatment arm
#' @param n.C [integer, >0] number of patients in the control arm
#' @param p1C [integer, >0] proportion of events of interest in the control group. Can be NULL if and only if \code{(b.1T, b.1C, b.2T, b.2C)}
#' are provided.
#' @param sHR [double, >0] pre-specified sub-distribution hazard ratio for event of interest
#' @param v.1C,v.1T,v.2C,v.2T shape parameters for distribution of time to event of interest in control/treatment (C/T) group and of time 
#' to competing event in control/treatment (C/T) group respectively
#' @param b.1C,b.1T,b.2C,b.2T rate parameters for distribution of time to event of interest in control/treatment (C/T) group and of time 
#' to competing event in control/treatment (C/T) group respectively. Can be NULL if and only if \code{(p.1C, sHR)} are provided.
#' @param cens.distrib [character] censoring distribution. Can be \code{"exponential"} for exponential censoring or \code{"uniform"} for
#' uniform censoring. NULL means no censoring.
#' @param param.cens [>0] parameter for censoring distribution. Should be a double for rate parameter of exponential censoring or a vector 
#' of doubles for lower and upper bounds of uniform censoring. NULL means no censoring
#' @param latent [logical] If \code{TRUE} also export the latent variables (e.g. true event times, event types and censoring times). NULL
#' sets this parameter to \code{FALSE}.
#' 
#' @details 
#' 
#' Arguments in the list \code{argsTTE}:
#'     \itemize{
#'     \item\code{CR} should competing risks be simulated? \cr 
#'     \item\code{rates.T} hazard corresponding to each endpoint (time to event endpoint, treatment group). \cr 
#'     \item\code{rates.C} same as \code{rates.T} but for the control group. \cr
#'     \item\code{rates.CR} same as \code{rates.T} but for the competing event (same in both groups). \cr
#'     \item\code{rates.Censoring} Censoring same as \code{rates.T} but for the censoring. \cr
#'     \item\code{name} names of the time to event variables. \cr
#'     \item\code{nameCensoring} names of the event type indicators. \cr
#'     }
#'     
#' @examples
#'
#' #### Providing p.1C and sHR ####
#' d <- simCompetingRisks(n.T = 5000, n.C = 5000, p.1C = 0.55, v.1C = -0.30, v.1T = -0.30, v.2C = -0.30, v.2T = -0.30, 
#' sHR = 0.5, b.1T = NULL, b.1C = NULL, b.2T = NULL, b.2C = NULL)
#' 
#' #### Providing the rate parameters ####
#' d <- simCompetingRisks(n.T=5000, n.C=5000, p.1C = NULL, v.1C = -0.30, v.1T = -0.30, v.2C = -0.30, v.2T = -0.30, 
#' sHR = NULL, b.1T = b.1T, b.1C = b.1C, b.2T = b.2T, b.2C = b.2C)
#' 
#' #### With exponential censoring ####
#' d <- simCompetingRisks(n.T = 5000, n.C = 5000, p.1C = 0.55, v.1C = -0.30, v.1T = -0.30, v.2C = -0.30, v.2T = -0.30, 
#' sHR = 0.5, b.1T = NULL, b.1C = NULL, b.2T = NULL, b.2C = NULL, cens.distrib = "exponential", param.cens = 0.8, 
#' latent = T)
#'
#' ### With uniform censoring ####
#' d <- simCompetingRisks(n.T = 5000, n.C = 5000, p.1C = 0.55, v.1C = -0.30, v.1T = -0.30, v.2C = -0.30, v.2T = -0.30, 
#' sHR = 0.5, b.1T = NULL, b.1C = NULL, b.2T = NULL, b.2C = NULL, cens.distrib = "uniform", param.cens = c(0, 7), 
#' latent=T)        
#' 
#' @keywords function simulations
#'

## * Function simCompetingRisks
#' @rdname simulation_CR
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
  u <- runif(n, 0, 1)
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
      data$censoring.time <- rexp(n, rate = param.cens[1])
    }
    else if (cens.distrib == "uniform") {
      if(is.na(param.cens[2])) {
        stop("Missing parameter for uniform censoring distribution")
      }
      data$censoring.time <- runif(n, min = param.cens[1], max = param.cens[2])
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
