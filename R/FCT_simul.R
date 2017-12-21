## * Documentation - simulBT
#' @name Simulation function
#' @rdname simulation
#' @title Simulation of data for the BuyseTest
#' @aliases simulBT
#' 
#' @description Simulate binary, continuous or time to event data, possibly with strata.
#' 
#' @param n.T number of patients in the treatment arm
#' @param n.C number of patients in the control arm
#' @param format the format of the output. Can be "data.table", "data.frame" or "matrix"
#' @param argsBin a list of arguments to be passed to simulBT_bin. They specify the distribution parameters of the binary endpoints
#' @param argsCont a list of arguments to be passed to simulBT_continous. They specify the distribution parameters of the continuous endpoints
#' @param argsTTE a list of arguments to be passed to simulBT_TTE. They specify the distribution parameters of the time to event endpoints
#' @param n.strata number of strata. \code{NULL} indicates no strata.
#' @param names.strata name of the strata varaibles. Must have same length as \code{n.strata}.
#' 
#' @details 
#' Built on the lvm and sim functions from the lava package.
#' 
#' Arguments in the list \code{argsBin}:
#' \itemize{
#'     \item\code{p.T} probability of event of each endpoint (binary endpoint, treatment group). \cr 
#'     \item\code{p.C} same as \code{p.T} but for the control group. \cr
#'     \item\code{name} names of the binary variables. \cr
#' }
#' 
#' Arguments in the list \code{argsCont}:
#'     \itemize{
#'     \item\code{mu.T} expected value of each endpoint (continuous endpoint, treatment group). \cr 
#'     \item\code{mu.C} same as \code{mu.C} but for the control group. \cr
#'     \item\code{sigma.T} standard deviation of the values of each endpoint (continuous endpoint, treatment group). \cr 
#'     \item\code{sigma.C} same as \code{sigma.T} but for the control group. \cr
#'     \item\code{name} names of the continous variables.
#'     }
#' 
#' Arguments in the list \code{argsTTE}:
#'     \itemize{
#'     \item\code{rates.T} hazard corresponding to each endpoint (time to event endpoint, treatment group). \cr 
#'     \item\code{rates.C} same as \code{rates.T} but for the control group. \cr
#'     \item\code{rates.Censoring} Censoring same as \code{rates.T} but for the censoring. \cr
#'     \item\code{name} names of the time to event variables. \cr
#'     \item\code{nameCensoring} names of the event type indicators. \cr
#'     }
#'     
#' @examples 
#' simulBT(1e3)
#' 
#' ## only binary endpoints
#' simulBT(1e3, argsBin = list(p.T = c(3:5/10)), argsCont = NULL, argsTTE = NULL)
#' 
#' ## only continous endpoints
#' simulBT(1e3, argsBin = NULL, argsCont = list(mu.T = c(3:5/10), sigma.T = rep(1,3)), argsTTE = NULL)
#' 
#' ## only TTE endpoints
#' simulBT(1e3, argsBin = NULL, argsCont = NULL, 
#'         argsTTE = list(rates.T = c(3:5/10), rates.Censoring = rep(1,3)))
#'         
#' ## one strata with 5 levels
#' simulBT(1e3, n.strata = 5)
#' simulBT(1e3, n.strata = 5, names.strata = "grade")
#' 
#' ## several strata
#' simulBT(1e3, n.strata = c(2,4), names.strata = c("Gender","AgeCategory"))
#' 
#' @keywords function simulations
#'

## * Function simulBT
#' @rdname simulation
#' @export
simulBT <- function(n.T, n.C = NULL, 
                    argsBin = list(), argsCont = list(), argsTTE = list(),
                    n.strata = NULL, names.strata = NULL, format = "data.table"){
  
  if(is.null(names.strata) && !is.null(n.strata)){
    if(length(n.strata)==1){names.strata <- "strata"}else{names.strata <- paste0("strataVar",1:n.strata)}
  }
  
  ## ** test
  if(is.null(n.C)){n.C <- n.T}
  validNumeric(n.C, min = 0, validLength = 1, method = "simulBT")
  validNumeric(n.T, min = 0, validLength = 1, method = "simulBT")
  validInteger(n.strata, validLength = NULL, refuse.NULL = FALSE, min = 1, method = "simulBT")
  validCharacter(format, validLength = 1, validValues = c("data.table","data.frame","matrix"), method = "simulBT")
  
  ## ** lvm 
  mT.lvm <- lvm()
  mC.lvm <- lvm()
  if(!is.null(argsBin)){
    newLVM <- do.call("simulBT_bin", args = c(list(modelT = mT.lvm, modelC = mC.lvm), argsBin))
    mT.lvm <- newLVM$modelT
    mC.lvm <- newLVM$modelC
  }
  if(!is.null(argsCont)){
    newLVM <- do.call("simulBT_cont", args = c(list(modelT = mT.lvm, modelC = mC.lvm), argsCont))
    mT.lvm <- newLVM$modelT
    mC.lvm <- newLVM$modelC
  }
  if(!is.null(argsTTE)){
    newLVM <- do.call("simulBT_TTE", args = c(list(modelT = mT.lvm, modelC = mC.lvm), argsTTE))
    mT.lvm <- newLVM$modelT
    mC.lvm <- newLVM$modelC
  }
  
  ## ** strata
  if(!is.null(n.strata)){
    validCharacter(names.strata, validLength = length(n.strata), refuse.NULL = TRUE, method = "simulBT")
    
    for(iterS in 1:length(n.strata)){
      if(any(names.strata[iterS] %in% lava::vars(mT.lvm))){
        stop("simulBT_bin: variable already in the LVM \n",
             "variable: ",paste(names.strata[iterS][names.strata[iterS] %in% lava::vars(mT.lvm)], collapse = " "),"\n")
      }
      
      lava::categorical(mT.lvm, labels = letters[1:n.strata[iterS]]) <- names.strata[iterS]
      lava::categorical(mC.lvm, labels = letters[1:n.strata[iterS]]) <- names.strata[iterS]   
    }
  }
  
  ## ** simulation
  df.T <- cbind(Treatment = 1, lava::sim(mT.lvm, n.T))
  df.C <- cbind(Treatment = 0, lava::sim(mC.lvm, n.C))
  if(!is.null(argsTTE)){
    df.T <- df.T[,setdiff(names(df.T),newLVM$latent)]
    df.C <- df.C[,setdiff(names(df.C),newLVM$latent)]
  }
  
  ## ** export
  res <- do.call(format, args =  rbind(df.T, df.C))
  return(res)
}

## * Function simulBT_bin
simulBT_bin <- function(modelT, modelC, p.T = 0.5, p.C = NULL, name = NULL){

  ## ** intialisation
  n.endpoints <- length(p.T)
  if(is.null(name)){ 
    if(n.endpoints == 1){name <- "toxicity"}else{name <- paste0("toxicity",1:n.endpoints)}
  }
  if(is.null(p.C)){p.C <- p.T}
  
  ## ** tests
  validNumeric(p.T, min = 0, max = 1, validLength = NULL, method = "simulBT")
  validNumeric(p.C, min = 0, max = 1, validLength = n.endpoints, method = "simulBT")
  validCharacter(name, validLength = n.endpoints, method = "simulBT")
  
  ## ** model
  for(iterE in 1:n.endpoints){
    if(any(name[iterE] %in% lava::vars(modelT))){
      stop("simulBT_bin: variable already in the LVM \n",
           "variable: ",paste(name[iterE][name[iterE] %in% lava::vars(modelT)], collapse = " "),"\n")
    }
    
    lava::distribution(modelT, name[iterE]) <- lava::binomial.lvm(link = "identity", p = p.T[iterE])
    lava::distribution(modelC, name[iterE]) <- lava::binomial.lvm(link = "identity", p = p.C[iterE])
  }
  
  ## ** export
  return(list(modelT = modelC, modelC = modelC))
  
}

## * Function simulBT_cont
simulBT_cont <- function(modelT, modelC, mu.T = 0, sigma.T = 1, mu.C = NULL, sigma.C = NULL, name = NULL){
  
  ## ** intialisation
  n.endpoints <- length(mu.T)
  if(is.null(name)){ 
    if(n.endpoints == 1){name <- "score"}else{name <- paste0("score",1:n.endpoints)}
  }
  if(is.null(mu.C)){mu.C <- mu.T}
  if(is.null(sigma.C)){sigma.C <- sigma.T}
  
  ## ** tests
  validNumeric(mu.T, validLength = NULL, method = "simulBT")
  validNumeric(sigma.T, validLength = n.endpoints, min = 0, method = "simulBT")
  validNumeric(mu.C, validLength = n.endpoints, method = "simulBT")
  validNumeric(sigma.C, validLength = n.endpoints, min = 0, method = "simulBT")
  validCharacter(name, validLength = n.endpoints, method = "simulBT")
  
  ## ** model
  for(iterE in 1:n.endpoints){
    if(any(name[iterE] %in% lava::vars(modelT))){
      stop("simulBT_cont: variable already in the LVM \n",
           "variable: ",paste(name[iterE][name[iterE] %in% lava::vars(modelT)], collapse = " "),"\n")
    }
    
    lava::distribution(modelT, name[iterE]) <- lava::gaussian.lvm(link = "identity", mean = mu.T[iterE], sd = sigma.T[iterE])
    lava::distribution(modelC, name[iterE]) <- lava::gaussian.lvm(link = "identity", mean = mu.C[iterE], sd = sigma.C[iterE])
  }
  
  ## ** export
  return(list(modelT = modelT, modelC = modelC))
}

## * Function simulBT_TTE
simulBT_TTE <- function(modelT, modelC, 
                        rates.T = 2, rates.C = NULL, rates.Censoring = 1, sigma.C = NULL, 
                        name = NULL, nameCensoring = NULL){
  
  ## ** intialisation
  n.endpoints <- length(rates.T)
  if(is.null(name)){ 
    if(n.endpoints == 1){name <- "eventtime"}else{name <- paste0("eventtime",1:n.endpoints)}
  }
  if(is.null(nameCensoring)){ 
    if(n.endpoints == 1){nameCensoring <- "status"}else{nameCensoring <- paste0("status",1:n.endpoints)}
  }
  if(is.null(rates.C)){rates.C <- rates.T}
  
  name0 <- paste0(name,"Uncensored")
  nameC <- paste0(name,"Censoring")
  
  ## ** tests
  validNumeric(rates.T, validLength = NULL, method = "simulBT")
  validNumeric(rates.C, validLength = n.endpoints, min = 0, method = "simulBT")
  validNumeric(rates.Censoring, validLength = n.endpoints, min = 0, method = "simulBT")
  validCharacter(name, validLength = n.endpoints, method = "simulBT")
  validCharacter(nameCensoring, validLength = n.endpoints, method = "simulBT")
  
  
  ## ** model
  for(iterE in 1:n.endpoints){
    allvarE <- c(name[iterE], name0[iterE], nameC[iterE], nameCensoring[iterE])
    if(any(allvarE %in% lava::vars(modelT))){
      stop("simulBT_TTE: variable already in the LVM \n",
           "variable: ",paste(allvarE[allvarE %in% lava::vars(modelT)], collapse = " "),"\n")
    }
    
    lava::distribution(modelT, name0[iterE]) <- lava::coxExponential.lvm(rate=rates.T[iterE])
    lava::distribution(modelT, nameC[iterE]) <- lava::coxExponential.lvm(rate=rates.Censoring[iterE])
    txtSurv <- paste0(name[iterE], "~min(",name0[iterE],"=1,",nameC[iterE],"=0)")
    modelT <- lava::eventTime(modelT, stats::as.formula(txtSurv), nameCensoring[iterE])
    
    lava::distribution(modelC, name0[iterE]) <- lava::coxExponential.lvm(rate=rates.C[iterE])
    lava::distribution(modelC, nameC[iterE]) <- lava::coxExponential.lvm(rate=rates.Censoring[iterE])
    txtSurv <- paste0(name[iterE], "~min(",name0[iterE],"=1,",nameC[iterE],"=0)")
    modelC <- lava::eventTime(modelC, stats::as.formula(txtSurv), nameCensoring[iterE])
  }
  
  ## ** export
  return(list(modelT = modelT, modelC = modelC, latent = c(name0,nameC)))
  
}
