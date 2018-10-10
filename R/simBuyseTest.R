## * Documentation - simBuyseTest
#' @name Simulation function
#' @rdname simulation
#' @title Simulation of data for the BuyseTest
#' @aliases simBuyseTest
#' 
#' @description Simulate binary, continuous or time to event data, possibly with strata.
#' Outcomes are simulated independently of each other and independently of the strata variable.
#' 
#' @param n.T [integer, >0] number of patients in the treatment arm
#' @param n.C [integer, >0] number of patients in the control arm
#' @param format [character] the format of the output. Can be \code{"data.table"}, \code{"data.frame"} or \code{"matrix"}.
#' @param argsBin [list] arguments to be passed to \code{simBuyseTest_bin}. They specify the distribution parameters of the binary endpoints.
#' @param argsCont [list] arguments to be passed to \code{simBuyseTest_continuous}. They specify the distribution parameters of the continuous endpoints.
#' @param argsTTE [list]  arguments to be passed to \code{simBuyseTest_TTE}. They specify the distribution parameters of the time to event endpoints.
#' @param n.strata [integer, >0] number of strata. \code{NULL} indicates no strata.
#' @param names.strata [character vector] name of the strata variables. Must have same length as \code{n.strata}.
#' @param latent [logical] If \code{TRUE} also export the latent variables (e.g. censoring times or event times).
#' 
#' @details 
#' This function is built upon the \code{lvm} and \code{sim} functions from the lava package.
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
#'     \item\code{name} names of the continuous variables.
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
#' n <- 1e2
#'
#' #### default option ####
#' simBuyseTest(n)
#'
#' ## with a strata variable having 5 levels
#' simBuyseTest(n, n.strata = 5)
#' ## with a strata variable named grade
#' simBuyseTest(n, n.strata = 5, names.strata = "grade")
#' ## several strata variables
#' simBuyseTest(1e3, n.strata = c(2,4), names.strata = c("Gender","AgeCategory"))
#' 
#' #### only binary endpoints ####
#' args <- list(p.T = c(3:5/10))
#' simBuyseTest(n, argsBin = args, argsCont = NULL, argsTTE = NULL)
#' 
#' #### only continuous endpoints ####
#' args <- list(mu.T = c(3:5/10), sigma.T = rep(1,3))
#' simBuyseTest(n, argsBin = NULL, argsCont = args, argsTTE = NULL)
#' 
#' #### only TTE endpoints ####
#' args <- list(rates.T = c(3:5/10), rates.Censoring = rep(1,3))
#' simBuyseTest(n, argsBin = NULL, argsCont = NULL, argsTTE = args)
#'         
#' 
#' @keywords function simulations
#'

## * Function simBuyseTest
#' @rdname simulation
#' @export
simBuyseTest <- function(n.T, n.C = NULL, 
                         argsBin = list(), argsCont = list(), argsTTE = list(),
                         n.strata = NULL, names.strata = NULL, format = "data.table",
                         latent = FALSE){

    option <- BuyseTest.options()
    if(is.null(names.strata) && !is.null(n.strata)){
        if(length(n.strata)==1){
            names.strata <- "strata"
        }else{
            names.strata <- paste0("strataVar",1:n.strata)
        }
    }
  
    ## ** check arguments
    if(is.null(n.C)){n.C <- n.T}

    if(option$check){
        validNumeric(n.C,
                     min = 0,
                     valid.length = 1,
                     method = "simBuyseTest")
        validNumeric(n.T,
                     min = 0,
                     valid.length = 1,
                     method = "simBuyseTest")
        validInteger(n.strata,
                     valid.length = NULL,
                     refuse.NULL = FALSE,
                     min = 1,
                     method = "simBuyseTest")
        validCharacter(format,
                       valid.length = 1,
                       valid.values = c("data.table","data.frame","matrix"),
                       method = "simBuyseTest")
    }
    
    ## ** build the generative model
    mT.lvm <- lvm()
    mC.lvm <- lvm()
    if(!is.null(argsBin)){
        newLVM <- do.call("simBuyseTest_bin", args = c(list(modelT = mT.lvm, modelC = mC.lvm, check = option$check), argsBin))
        mT.lvm <- newLVM$modelT
        mC.lvm <- newLVM$modelC
    }
    if(!is.null(argsCont)){
        newLVM <- do.call("simBuyseTest_cont", args = c(list(modelT = mT.lvm, modelC = mC.lvm, check = option$check), argsCont))
        mT.lvm <- newLVM$modelT
        mC.lvm <- newLVM$modelC
    }
    if(!is.null(argsTTE)){
        newLVM <- do.call("simBuyseTest_TTE", args = c(list(modelT = mT.lvm, modelC = mC.lvm, check = option$check), argsTTE))
        mT.lvm <- newLVM$modelT
        mC.lvm <- newLVM$modelC
    }
  
    ## ** add strata variable to the generative model
    if(!is.null(n.strata)){
        if(option$check){
            validCharacter(names.strata,
                           valid.length = length(n.strata),
                           refuse.NULL = TRUE,
                           method = "simBuyseTest")
        }
        
        for(iterS in 1:length(n.strata)){
            if(any(names.strata[iterS] %in% lava::vars(mT.lvm))){
                stop("simBuyseTest: variable already in the LVM \n",
                     "variable: ",paste(names.strata[iterS][names.strata[iterS] %in% lava::vars(mT.lvm)], collapse = " "),"\n")
            }
      
            lava::categorical(mT.lvm, labels = letters[1:n.strata[iterS]]) <- names.strata[iterS]
            lava::categorical(mC.lvm, labels = letters[1:n.strata[iterS]]) <- names.strata[iterS]   
        }
    }
  
    ## ** simulate data from the generative model
    df.T <- cbind(Treatment = "T", lava::sim(mT.lvm, n.T, latent = latent))
    df.C <- cbind(Treatment = "C", lava::sim(mC.lvm, n.C, latent = latent))
  
  
    ## ** export
    res <- do.call(format, args =  rbind(df.C, df.T))
    return(res)
}

## * Function simBuyseTest_bin
simBuyseTest_bin <- function(modelT, modelC, p.T = 0.5, p.C = NULL, name = NULL, check){

  ## ** initialisation
  n.endpoints <- length(p.T)
  if(is.null(name)){ 
    if(n.endpoints == 1){name <- "toxicity"}else{name <- paste0("toxicity",1:n.endpoints)}
  }
  if(is.null(p.C)){p.C <- p.T}
  
    ## ** tests
    if(check){
        validNumeric(p.T,
                     min = 0,
                     max = 1,
                     valid.length = NULL,
                     method = "simBuyseTest")
        validNumeric(p.C,
                     min = 0,
                     max = 1,
                     valid.length = n.endpoints,
                     method = "simBuyseTest")
        validCharacter(name,
                       valid.length = n.endpoints,
                       method = "simBuyseTest")
    }
    
    ## ** model
    for(iterE in 1:n.endpoints){
        if(any(name[iterE] %in% lava::vars(modelT))){
            stop("simBuyseTest_bin: variable already in the LVM \n",
                 "variable: ",paste(name[iterE][name[iterE] %in% lava::vars(modelT)], collapse = " "),"\n")
        }
    
    lava::distribution(modelT, name[iterE]) <- lava::binomial.lvm(link = "identity", p = p.T[iterE])
    lava::distribution(modelC, name[iterE]) <- lava::binomial.lvm(link = "identity", p = p.C[iterE])
  }
  
  ## ** export
  return(list(modelT = modelC, modelC = modelC))
  
}

## * Function simBuyseTest_cont
simBuyseTest_cont <- function(modelT,
                              modelC,
                              mu.T = 0,
                              sigma.T = 1,
                              mu.C = NULL,
                              sigma.C = NULL,
                              name = NULL,
                              check){
    
    ## ** initialisation
    n.endpoints <- length(mu.T)
    if(is.null(name)){ 
        if(n.endpoints == 1){name <- "score"}else{name <- paste0("score",1:n.endpoints)}
    }
    if(is.null(mu.C)){mu.C <- mu.T}
    if(is.null(sigma.C)){sigma.C <- sigma.T}
    
    ## ** tests
    if(check){
        validNumeric(mu.T,
                     valid.length = NULL,
                     method = "simBuyseTest")
        validNumeric(sigma.T,
                     valid.length = n.endpoints,
                     min = 0,
                     method = "simBuyseTest")
        validNumeric(mu.C,
                     valid.length = n.endpoints,
                     method = "simBuyseTest")
        validNumeric(sigma.C,
                     valid.length = n.endpoints,
                     min = 0,
                     method = "simBuyseTest")
        validCharacter(name,
                       valid.length = n.endpoints,
                       method = "simBuyseTest")
    }
    
    ## ** model
    for(iterE in 1:n.endpoints){
        if(any(name[iterE] %in% lava::vars(modelT))){
            stop("simBuyseTest_cont: variable already in the LVM \n",
                 "variable: ",paste(name[iterE][name[iterE] %in% lava::vars(modelT)], collapse = " "),"\n")
        }
        
        lava::distribution(modelT, name[iterE]) <- lava::gaussian.lvm(link = "identity",
                                                                      mean = mu.T[iterE],
                                                                      sd = sigma.T[iterE])
        lava::distribution(modelC, name[iterE]) <- lava::gaussian.lvm(link = "identity",
                                                                      mean = mu.C[iterE],
                                                                      sd = sigma.C[iterE])
    }
    
    ## ** export
    return(list(modelT = modelT, modelC = modelC))
}

## * Function simBuyseTest_TTE
simBuyseTest_TTE <- function(modelT,
                             modelC, 
                             rates.T = 2,
                             rates.C = NULL,
                             rates.Censoring = 1,
                             sigma.C = NULL, 
                             name = NULL,
                             nameCensoring = NULL,
                             check){
    
    ## ** initialisation
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
    if(check){
        validNumeric(rates.T,
                     valid.length = NULL,
                     method = "simBuyseTest")
        validNumeric(rates.C,
                     valid.length = n.endpoints,
                     min = 0,
                     method = "simBuyseTest")
        validNumeric(rates.Censoring,
                     valid.length = n.endpoints,
                     min = 0,
                     method = "simBuyseTest")
        validCharacter(name,
                       valid.length = n.endpoints,
                       method = "simBuyseTest")
        validCharacter(nameCensoring,
                       valid.length = n.endpoints,
                       method = "simBuyseTest")  
    }
    
    ## ** model
    for(iterE in 1:n.endpoints){
        allvarE <- c(name[iterE], name0[iterE], nameC[iterE], nameCensoring[iterE])
        if(any(allvarE %in% lava::vars(modelT))){
            stop("simBuyseTest_TTE: variable already in the LVM \n",
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

        formula.latent <- as.formula(paste0("~",name0[iterE],"+",nameC[iterE]))
        latent(modelT) <- formula.latent
        latent(modelC) <- formula.latent
    }

    ## ** export
    return(list(modelT = modelT, modelC = modelC, latent = c(name0,nameC)))
    
}
