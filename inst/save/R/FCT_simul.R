#' @name Simulation function
#' @rdname simulation
#' @title Simulation of data for the BuyseTest
#' @aliases simulBT
#' 
#' @description Simulate binary, continuous or time to event data, possibly with strata.
#' 
#' @param n.T number of patients in the treatment arm
#' @param n.C number of patients in the control arm
#' @param n.strata number of strata. \code{NULL} indicates no strata
#' @param format the format of the output. Can be "data.table", "data.frame" or "matrix"
#' @param argsBin a list of arguments to be passed to simulBT_bin. They specify the distribution parameters of the binary endpoints
#' @param argsCont a list of arguments to be passed to simulBT_continous. They specify the distribution parameters of the continuous endpoints
#' @param argsTTE a list of arguments to be passed to simulBT_TTE. They specify the distribution parameters of the time to event endpoints
#' @param p.T probability of event of each endpoint (binary endpoint, treatment group)
#' @param p.C same as \code{p.T} but for the control group
#' @param mu.T expected value of each endpoint (continuous endpoint, treatment group)
#' @param mu.C same as \code{mu.C} but for the control group
#' @param sigma.T standard deviation of the values of each endpoint (continuous endpoint, treatment group)
#' @param sigma.C same as \code{sigma.T} but for the control group
#' @param rates.T hazard corresponding to each endpoint (time to event endpoint, treatment group)
#' @param rates.C same as \code{rates.T} but for the control group
#' @param rates.Censor same as \code{rates.T} but for the censoring.
#' 
#' @details 
#' Built on the lvm and sim functions from the lava package.
#' 
#' @examples 
#' simulBT(1e3)
#' simulBT(1e3, argsBin = list(p.T = c(3:5/10)), argsCont = NULL, argsTTE = NULL)
#' simulBT(1e3, argsBin = NULL, argsCont = list(mu.T = c(3:5/10), sigma.T = rep(1,3)), argsTTE = NULL)
#' simulBT(1e3, argsBin = NULL, argsCont = NULL, argsTTE = list(rates.T = c(3:5/10), rates.Censor = rep(1,3)))
#' 
#' @keywords function simulations
#' 

#' @rdname simulation
#' @export
simulBT <- function(n.T, n.C = NULL, n.strata = NULL, format = "data.table",
                    argsBin = list(), argsCont = list(), argsTTE = list()){
  
  res <- NULL
  if(!is.null(argsBin)){
    resBin <- do.call("simulBT_bin", args = c(list(n.T = n.T, n.C = n.C, n.strata = n.strata, format = "data.table"),
                                    argsBin))
    res <- resBin
  }
  if(!is.null(argsCont)){
    resCont <- do.call("simulBT_cont", args = c(list(n.T = n.T, n.C = n.C, n.strata = n.strata, format = "data.table"),
                                    argsCont))
    if(is.null(res)){res <- resCont}else{res <- cbind(res, resCont[,"Treatment" := NULL])}
  }
  if(!is.null(argsTTE)){
    resTTE <- do.call("simulBT_TTE", args = c(list(n.T = n.T, n.C = n.C, n.strata = n.strata, format = "data.table"),
                                    argsTTE))
    if(is.null(res)){res <- resTTE}else{res <- cbind(res, resTTE[,"Treatment" := NULL])}
  }
  
  res <- do.call(format, args = res)
  return(res)
}

#' @rdname simulation
#' @export
simulBT_bin <- function(n.T, n.C = NULL, p.T = 0.5, p.C = NULL, n.strata = NULL, format =  "data.table"){

  
  if(is.null(n.C)){n.C <- n.T}
  if(is.null(p.C)){p.C <- p.T}
  
  validNumeric(n.C, min = 0, validLength = 1, method = "simulBT")
  validNumeric(n.T, min = 0, validLength = 1, method = "simulBT")
  validCharacter(format, validLength = 1, validValues = c("data.table","data.frame","matrix"), method = "simulBT")
  
  validNumeric(p.T, min = 0, max = 1, validLength = NULL, method = "simulBT")
  n.endpoints <- length(p.T)
  validNumeric(p.C, min = 0, max = 1, validLength = n.endpoints, method = "simulBT")
  validNumeric(n.strata, validLength = 1, refuse.NULL = FALSE, min = 1, method = "simulBT")
  
  #### model
  formula.lvm <- lapply(paste0(paste0("Y_bin",1:n.endpoints),"~1", sep = ""),as.formula)
  simul <- lava::lvm(formula.lvm)
    
  ## Treatment
  simul.T <- simul
  for(iterE in 1:n.endpoints){
    lava::distribution(simul.T, paste0("Y_bin",iterE)) <- lava::binomial.lvm(link = "identity", p = p.T[iterE])
  }
  if(!is.null(n.strata)){
    lava::categorical(simul.T, labels = letters[1:n.strata]) <- "strata" 
  }
  df.T <- cbind(Treatment = 1, lava::sim(simul.T, n.T))
  
  ## Control
  simul.C <- simul
  for(iterE in 1:n.endpoints){
    lava::distribution(simul.C, paste0("Y_bin",iterE)) <- lava::binomial.lvm(link = "identity", p = p.C[iterE])
  }
  if(!is.null(n.strata)){
    lava::categorical(simul.C, labels = letters[1:n.strata]) <- "strata" 
  }
  df.C <- cbind(Treatment = 0, lava::sim(simul.C, n.C))
  
  ## Merge
  df <- do.call(format, args = rbind(df.T, df.C))
  return(df)
 
}

#' @rdname simulation
#' @export
simulBT_cont <- function(n.T, n.C = NULL, mu.T = 0, sigma.T = 1, mu.C = NULL, sigma.C = NULL, n.strata = NULL, format =  "data.table"){
  
  if(is.null(n.C)){n.C <- n.T}
  if(is.null(mu.C)){mu.C <- mu.T}
  if(is.null(sigma.C)){sigma.C <- sigma.T}
  
  validNumeric(n.C, min = 0, validLength = 1, method = "simulBT")
  validNumeric(n.T, min = 0, validLength = 1, method = "simulBT")
  validCharacter(format, validLength = 1, validValues = c("data.table","data.frame","matrix"), method = "simulBT")
  
  validNumeric(mu.T, validLength = NULL, method = "simulBT")
  n.endpoints <- length(mu.T)
  validNumeric(sigma.T, validLength = n.endpoints, min = 0, method = "simulBT")
  validNumeric(mu.C, validLength = n.endpoints, method = "simulBT")
  validNumeric(sigma.C, validLength = n.endpoints, min = 0, method = "simulBT")
  validNumeric(n.strata, validLength = 1, refuse.NULL = FALSE, min = 1, method = "simulBT")
  
  #### model
  formula.lvm <- lapply(paste0(paste0("Y_cont",1:n.endpoints),"~1", sep = ""),as.formula)
  simul <- lava::lvm(formula.lvm)
  
  ## Treatment
  simul.T <- simul
  for(iterE in 1:n.endpoints){
    lava::distribution(simul.T, paste0("Y_cont",iterE)) <- lava::gaussian.lvm(link = "identity", mean = mu.T[iterE], sd = sigma.T[iterE])
  }
  if(!is.null(n.strata)){
    lava::categorical(simul.T, labels = letters[1:n.strata]) <- "strata" 
  }
  df.T <- cbind(Treatment = 1, lava::sim(simul.T, n.T))
  
  ## Control
  simul.C <- simul
  for(iterE in 1:n.endpoints){
    lava::distribution(simul.C, paste0("Y_cont",iterE)) <- lava::gaussian.lvm(link = "identity", mean = mu.C[iterE], sd = sigma.C[iterE])
  }
  if(!is.null(n.strata)){
    lava::categorical(simul.C, labels = letters[1:n.strata]) <- "strata" 
  }
  df.C <- cbind(Treatment = 0, lava::sim(simul.C, n.C))
  
  ## Merge
  df <- do.call(format, args = rbind(df.T, df.C))
  return(df)
  
}

#' @rdname simulation
#' @export
simulBT_TTE <- function(n.T, n.C = NULL, rates.T = 2, rates.C = NULL, rates.Censor = 1, n.strata = NULL, sigma.C = NULL, format =  "data.table"){
  
  if(is.null(n.C)){n.C <- n.T}
  if(is.null(rates.C)){rates.C <- rates.T}
  
  validNumeric(n.C, min = 0, validLength = 1, method = "simulBT")
  validNumeric(n.T, min = 0, validLength = 1, method = "simulBT")
  validCharacter(format, validLength = 1, validValues = c("data.table","data.frame","matrix"), method = "simulBT")
  
  validNumeric(rates.T, validLength = NULL, method = "simulBT")
  n.endpoints <- length(rates.T)
  validNumeric(rates.C, validLength = n.endpoints, min = 0, method = "simulBT")
  validNumeric(rates.Censor, validLength = n.endpoints, min = 0, method = "simulBT")
  validNumeric(n.strata, validLength = 1, refuse.NULL = FALSE, min = 1, method = "simulBT")
  
  #### model
  formula.lvm <- lapply(paste0(paste0("Y_TTE",1:n.endpoints),"~1", sep = ""),as.formula)
  simul <- lava::lvm(formula.lvm)
  
  ## Treatment
  simul.T <- simul
  
  for(iterE in 1:n.endpoints){
    lava::distribution(simul.T, paste0("T_TTE",iterE)) <- lava::coxExponential.lvm(rate=rates.T[iterE])
    lava::distribution(simul.T, paste0("Cens_TTE",iterE)) <- lava::coxExponential.lvm(rate=rates.Censor[iterE])
    txtSurv <- paste0(paste0("Y_TTE",iterE), "~min(",paste0("T_TTE",iterE),"=1,",paste0("Cens_TTE",iterE),"=0)")
    simul.T <- lava::eventTime(simul.T, as.formula(txtSurv), paste0("event",iterE))
  }
  if(!is.null(n.strata)){
    lava::categorical(simul.T, labels = letters[1:n.strata]) <- "strata" 
  }
  df.T <- cbind(Treatment = 1, lava::sim(simul.T, n.T))
  
  ## Control
  simul.C <- simul
  for(iterE in 1:n.endpoints){
    lava::distribution(simul.C, paste0("T_TTE",iterE)) <- lava::coxExponential.lvm(rate=rates.C[iterE])
    lava::distribution(simul.C, paste0("Cens_TTE",iterE)) <- lava::coxExponential.lvm(rate=rates.Censor[iterE])
    txtSurv <- paste0(paste0("Y_TTE",iterE), "~min(",paste0("T_TTE",iterE),"=1,",paste0("Cens_TTE",iterE),"=0)")
    simul.C <- lava::eventTime(simul.C, as.formula(txtSurv), paste0("event",iterE))
  }
  if(!is.null(n.strata)){
    lava::categorical(simul.C, labels = letters[1:n.strata]) <- "strata" 
  }
  df.C <- cbind(Treatment = 0, lava::sim(simul.C, n.C))
  
  ## Merge
  df <- do.call(format, args = rbind(df.T, df.C))
  return(df)
  
}
