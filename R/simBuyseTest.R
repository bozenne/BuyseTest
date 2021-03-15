## * Documentation - simBuyseTest
#' @name Simulate endpoints for GPC
#' @rdname simulation
#' @title Simulation of data for the BuyseTest
#' 
#' @description Simulate categorical, continuous or time to event endpoints, possibly along with a strata variable.
#' Categorical endpoints are simulated by thresholding a latent Gaussian variable (tobit model),
#' continuous endpoints are simulated using a Gaussian distribution,
#' and time to event endpoints are simulated using Weibull distributions for the event of interest, competing events, and censoring.
#' This function is built upon the \code{lvm} and \code{sim} functions from the lava package.
#' 
#' @param n.T [integer, >0] number of patients in the treatment arm
#' @param n.C [integer, >0] number of patients in the control arm
#' @param format [character] the format of the output. Can be \code{"data.table"}, \code{"data.frame"} or \code{"matrix"}.
#' @param argsBin [list] arguments to be passed to \code{simBuyseTest_bin}. They specify the distribution parameters of the categorical endpoints.
#' @param argsCont [list] arguments to be passed to \code{simBuyseTest_continuous}. They specify the distribution parameters of the continuous endpoints.
#' @param argsTTE [list]  arguments to be passed to \code{simBuyseTest_TTE}. They specify the distribution parameters of the time to event endpoints.
#' @param n.strata [integer, >0] number of strata. \code{NULL} indicates no strata.
#' @param names.strata [character vector] name of the strata variables. Must have same length as \code{n.strata}.
#' @param latent [logical] If \code{TRUE} also export the latent variables (e.g. censoring times or event times).
#' 
#' @details 
#' Endpoints are simulated independently of the strata variable and independently of each other,
#' with the exception of categorical endpoint and the time to event endpoints that can be correlated
#' by specifying a non-0 value for the \code{rho.T} and \code{rho.C} elements of the argument \code{argsBin}.
#' 
#' Arguments in the list \code{argsBin}:
#' \itemize{
#'     \item\code{p.T} list of probabilities for the values taken by each endpoint (categorical endpoint, treatment group). 
#'     \item\code{p.C} same as \code{p.T} but for the control group. 
#'     \item\code{rho.T} value of the regression coefficient between the underlying latent variable and the survival time. 
#'     \item\code{rho.C} same as \code{rho.T} but for the control group. 
#'     \item\code{name} names of the binary variables.
#' }
#' 
#' Arguments in the list \code{argsCont}:
#'     \itemize{
#'     \item\code{mu.T} expected value of each endpoint (continuous endpoint, treatment group). 
#'     \item\code{mu.C} same as \code{mu.C} but for the control group. 
#'     \item\code{sigma.T} standard deviation of the values of each endpoint (continuous endpoint, treatment group). 
#'     \item\code{sigma.C} same as \code{sigma.T} but for the control group. 
#'     \item\code{name} names of the continuous variables.
#'     }
#' 
#' Arguments in the list \code{argsTTE}:
#'     \itemize{
#'     \item\code{CR} should competing risks be simulated? 
#'     \item\code{scale.T,scale.C,scale.CR,scale.Censoring.T,scale.Censoring.C} scale parameter of the Weibull distribution for, respectively,
#'      the event of interest in the treatment group,
#'      the event of interest in the control group,
#'      the competing event in both groups,
#'      the censoring mechanism in the treatment group,
#'      the censoring mechanism in the control group
#'     \item\code{shape.T,shape.C,shape.CR,shape.Censoring.T,shape.Censoring.C} shape parameter of the Weibull distribution for, respectively,
#'      the event of interest in the treatment group,
#'      the event of interest in the control group,
#'      the competing event in both groups,
#'      the censoring mechanism in the treatment group,
#'      the censoring mechanism in the control group
#'     \item\code{name} names of the time to event variables. 
#'     \item\code{nameCensoring} names of the event type indicators. 
#'     }
#'     
#' @examples
#' library(data.table)
#' 
#' n <- 1e2
#'
#' #### by default ####
#' simBuyseTest(n)
#'
#' ## with a strata variable having 5 levels
#' simBuyseTest(n, n.strata = 5)
#' ## with a strata variable named grade
#' simBuyseTest(n, n.strata = 5, names.strata = "grade")
#' ## several strata variables
#' simBuyseTest(1e3, n.strata = c(2,4), names.strata = c("Gender","AgeCategory"))
#' 
#' #### only categorical endpoints ####
#' args <- list(p.T = list(c(low=0.1,moderate=0.5,high=0.4)))
#' dt.bin <- simBuyseTest(n, argsBin = args, argsCont = NULL, argsTTE = NULL)
#' table(dt.bin$toxicity)/NROW(dt.bin)
#' 
#' args <- list(p.T = list(c(low=0.1,moderate=0.5,high=0.4), c(0.1,0.9)))
#' dt.bin <- simBuyseTest(n, argsBin = args, argsCont = NULL, argsTTE = NULL)
#' table(dt.bin$toxicity1)/NROW(dt.bin)
#' table(dt.bin$toxicity2)/NROW(dt.bin)
#' 
#' #### only continuous endpoints ####
#' args <- list(mu.T = c(3:5/10), sigma.T = rep(1,3))
#' dt.cont <- simBuyseTest(n, argsBin = NULL, argsCont = args, argsTTE = NULL)
#' c(mean(dt.cont$score1), mean(dt.cont$score2), mean(dt.cont$score3))
#' c(sd(dt.cont$score1), sd(dt.cont$score2), sd(dt.cont$score3))
#' 
#' #### only TTE endpoints ####
#' args <- list(scale.T = c(3:5/10), scale.Censoring.T = rep(1,3))
#' dt.tte <- simBuyseTest(n, argsBin = NULL, argsCont = NULL, argsTTE = args)
#' 1/c(sum(dt.tte$eventtime1)/sum(dt.tte$status1),
#'   sum(dt.tte$eventtime2)/sum(dt.tte$status2),
#'   sum(dt.tte$eventtime3)/sum(dt.tte$status3))
#'         
#' 1/c(sum(dt.tte$eventtime1)/sum(dt.tte$status1==0),
#'   sum(dt.tte$eventtime2)/sum(dt.tte$status2==0),
#'   sum(dt.tte$eventtime3)/sum(dt.tte$status3==0))
#'
#' #### correlated categorical / time to event endpoint ####
#' args.bin <- list(p.T = list(c(low=0.1,moderate=0.5,high=0.4)), rho.T = 1)
#' args.tte <- list(scale.T = 2, scale.Censoring.T = 1)
#' dt.corr <- simBuyseTest(n, argsBin = args.bin, argsCont = NULL, argsTTE = args.tte)
#' 
#' 1/(sum(dt.corr$eventtime)/sum(dt.corr$status))
#' 1/(sum(dt.corr$eventtime)/sum(dt.corr$status==0))
#' table(dt.corr$toxicity)/NROW(dt.corr)
#' 
#' boxplot(eventtime ~ toxicity, data = dt.corr)
#' 
#' @keywords function simulations
#' @author Brice Ozenne

## * Function simBuyseTest
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
    lava::categorical(mT.lvm,labels=c("T")) <- "treatment"
    lava::categorical(mC.lvm,labels=c("C")) <- "treatment"
    if(!is.null(argsTTE)){
        newLVM <- do.call("simBuyseTest_TTE", args = c(list(modelT = mT.lvm, modelC = mC.lvm, check = option$check), argsTTE))
        mT.lvm <- newLVM$modelT
        mC.lvm <- newLVM$modelC
        latentTTE <- newLVM$latent0
        scale.T <- newLVM$scale.T
        scale.C <- newLVM$scale.C
        shape.T <- newLVM$shape.T
        shape.C <- newLVM$shape.C
    }else{
        latentTTE <- NULL
        scale.T <- NULL
        scale.C <- NULL
        shape.T <- NULL
        shape.C <- NULL
    }
    if(!is.null(argsBin)){
        newLVM <- do.call("simBuyseTest_bin", args = c(list(modelT = mT.lvm, modelC = mC.lvm, check = option$check,
                                                            latentTTE = latentTTE,
                                                            scale.T = scale.T, scale.C = scale.C, shape.T = shape.T, shape.C = shape.C),
                                                       argsBin))
        mT.lvm <- newLVM$modelT
        mC.lvm <- newLVM$modelC
    }
    if(!is.null(argsCont)){
        newLVM <- do.call("simBuyseTest_cont", args = c(list(modelT = mT.lvm, modelC = mC.lvm,check = option$check),
                                                        argsCont))
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
    df.T <- lava::sim(mT.lvm, n.T, latent = latent)
    df.C <- lava::sim(mC.lvm, n.C, latent = latent)
  
    ## ** export
    res <- do.call(format, args =  rbind(df.C, df.T))
    return(res)
}

## * Function simBuyseTest_bin
simBuyseTest_bin <- function(modelT,
                             modelC,
                             check,
                             latentTTE,
                             scale.T,
                             scale.C,
                             shape.T,
                             shape.C,
                             p.T = c("yes" = 0.5, "no" = 0.5),
                             p.C = NULL,
                             rho.T = NULL,
                             rho.C = NULL,
                             name = NULL){

    ## ** initialisation
    if(!is.null(p.T) && !is.list(p.T)){
        p.T <- list(p.T)
    }
    n.endpoints <- length(p.T)
    if(is.null(name)){ 
        if(n.endpoints == 1){name <- "toxicity"}else{name <- paste0("toxicity",1:n.endpoints)}
    }
    if(is.null(p.C)){
        p.C <- p.T
    }else if(!is.list(p.C)){
        p.C <- list(p.C)
    }
    if(is.null(rho.T)){
        rho.T <- rep(0, n.endpoints)
    }
    if(is.null(rho.C)){
        rho.C <- rho.T
    }
    names.values <- vector(mode = "list", length = n.endpoints)
    for(iterE in 1:n.endpoints){
        if(is.null(names(p.T[[iterE]]))){
            names.values[[iterE]] <- 1:length(p.T[[iterE]])
        }else{
            names.values[[iterE]] <- names(p.T[[iterE]])
        }
    }

    ## ** tests
    if(check){
        if(length(p.C)!=length(p.T)){
            stop("Arguments \'p.C\' and \'p.T\' must be a list with the same number of elements. \n",
                 "(each element defines the probability distribution of an endpoint; there must be the same number of endpoints in both groups) \n")
        }
        if(n.endpoints!=length(name)){
            stop("The length of arguments \'name\' does not match the number of endpoints defined by argument \'p.T\' \n")
        }
        validNumeric(rho.T,
                     valid.length = n.endpoints,
                     method = "simBuyseTest")
        validNumeric(rho.C,
                     valid.length = n.endpoints,
                     method = "simBuyseTest")
        if((any(rho.T!=0) || any(rho.C!=0)) && (n.endpoints != length(latentTTE))){
            stop("The number of time to event endpoints must match the number of categorical endpoints. \n")
        }
        for(iterE in 1:n.endpoints){
            validNumeric(p.T[[iterE]],
                         min = 0,
                         max = 1,
                         valid.length = NULL,
                         method = "simBuyseTest")
            if(sum(p.T[[iterE]])!=1){
                stop("For each endpoint, the sum of the probabilities in argument \'p.T\' must be 1. \n")
            }
            validNumeric(p.C[[iterE]],
                         min = 0,
                         max = 1,
                         valid.length = length(p.T[[iterE]]),
                         method = "simBuyseTest")
            if(sum(p.C[[iterE]])!=1){
                stop("For each endpoint, the sum of the probabilities in argument \'p.C\' must be 1. \n")
            }

            if(!identical(names(p.T[[iterE]]),names(p.C[[iterE]]))){
                stop("The names in arguments \'p.T\' and \'p.C\' must be the same. \n")
            }
        }
    }
    
    ## ** model
    for(iterE in 1:n.endpoints){
        if(any(name[iterE] %in% lava::vars(modelT))){
            stop("simBuyseTest_bin: variable already in the LVM \n",
                 "variable: ",paste(name[iterE][name[iterE] %in% lava::vars(modelT)], collapse = " "),"\n")
        }
        iLatent.T <- paste0("eta_",name[iterE])
        iLatent.C <- paste0("eta_",name[iterE])
        iCut.T <- qnormweibull(cumsum(p.T[[iterE]])[-length(p.T[[iterE]])], scale = scale.T[iterE], shape = shape.T[iterE], rho = rho.T[iterE])
        iFct.T <- paste0("function(x, xcut = c(",paste0(iCut.T,collapse=","),"), xname = c(\"",paste0(names.values[[iterE]],collapse="\",\""),"\")){\n",
                         "    return(factor(findInterval(x[,1], vec = xcut), levels = 0:length(xcut), labels = xname))\n",
                         "}")
        if(abs(rho.T[iterE]) > 1e-12){
            lava::regression(modelT) <- as.formula(paste0(iLatent.T," ~ ",rho.T," * ",latentTTE[iterE]))
        }
        modelT <- lava::`transform<-`(modelT, as.formula(paste0(name[iterE],"~",iLatent.T)), value = eval(parse(text = iFct.T)))
        lava::latent(modelT) <- as.formula(paste0("~",iLatent.T))

        iCut.C <- qnormweibull(cumsum(p.C[[iterE]])[-length(p.C[[iterE]])], scale = scale.C[iterE], shape = shape.C[iterE], rho = rho.C[iterE])
        iFct.C <- paste0("function(x, xcut = c(",paste0(iCut.C,collapse=","),"), xname = c(\"",paste0(names.values[[iterE]],collapse="\",\""),"\")){\n",
                       "    return(factor(findInterval(x[,1], vec = xcut), levels = 0:length(xcut), labels = xname))\n",
                       "}")
        if(abs(rho.C[iterE]) > 1e-12){
            lava::regression(modelC) <- as.formula(paste0(iLatent.C," ~ ",rho.C," * ",latentTTE[iterE]))
        }
        modelC <- lava::`transform<-`(modelC, as.formula(paste0(name[iterE],"~",iLatent.C)), value = eval(parse(text = iFct.C)))
        lava::latent(modelC) <- as.formula(paste0("~",iLatent.C))
        
    }
    
    ## ** export
    return(list(modelT = modelT, modelC = modelC))
    
}

## * Function simBuyseTest_cont
simBuyseTest_cont <- function(modelT,
                              modelC,
                              check,
                              mu.T = 0,
                              sigma.T = 1,
                              mu.C = NULL,
                              sigma.C = NULL,
                              name = NULL){
    
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
                             CR = FALSE,
                             scale.T = 1/2,
                             scale.C = NULL,
                             scale.CR = NULL,
                             shape.T = rep(1, length(scale.T)),
                             shape.C = NULL,
                             shape.CR = NULL,
                             scale.Censoring.T = rep(1, length(scale.T)),
                             scale.Censoring.C = NULL,
                             shape.Censoring.T = rep(1, length(scale.T)),
                             shape.Censoring.C = NULL,
                             name = NULL,
                             nameCensoring = NULL,
                             check){
    
    ## ** initialisation
    n.endpoints <- length(scale.T)
    if(is.null(name)){ 
        if(n.endpoints == 1){name <- "eventtime"}else{name <- paste0("eventtime",1:n.endpoints)}
    }
    if(is.null(nameCensoring)){ 
        if(n.endpoints == 1){nameCensoring <- "status"}else{nameCensoring <- paste0("status",1:n.endpoints)}
    }
    if(is.null(scale.C)){scale.C <- scale.T}
    if(is.null(scale.CR)){scale.CR <- scale.T}
    if(is.null(scale.Censoring.C)){scale.Censoring.C <- scale.Censoring.T}

    if(is.null(shape.C)){shape.C <- shape.T}
    if(is.null(shape.CR)){shape.CR <- shape.T}
    if(is.null(shape.Censoring.C)){shape.Censoring.C <- shape.Censoring.T}

    name0 <- paste0(name,"Uncensored")
    if(CR){
        nameCR <- paste0(name,"CompetingRisk")
    }
    nameC <- paste0(name,"Censoring")

    ## ** tests
    if(check){
        validLogical(CR,
                     valid.length = 1,
                     method = "simBuyseTest")
        
        validNumeric(scale.T,
                     valid.length = NULL,
                     method = "simBuyseTest")
        validNumeric(scale.C,
                     valid.length = n.endpoints,
                     min = 0,
                     method = "simBuyseTest")
        if(CR){
            validNumeric(scale.CR,
                         valid.length = n.endpoints,
                         min = 0,
                         method = "simBuyseTest")
        }
        validNumeric(scale.Censoring.T,
                     valid.length = n.endpoints,
                     min = 0,
                     method = "simBuyseTest")
        validNumeric(scale.Censoring.C,
                     valid.length = n.endpoints,
                     min = 0,
                     method = "simBuyseTest")

        validNumeric(shape.T,
                     valid.length = NULL,
                     method = "simBuyseTest")
        validNumeric(shape.C,
                     valid.length = n.endpoints,
                     min = 0,
                     method = "simBuyseTest")
        if(CR){
            validNumeric(shape.CR,
                         valid.length = n.endpoints,
                         min = 0,
                         method = "simBuyseTest")
        }
        validNumeric(shape.Censoring.T,
                     valid.length = n.endpoints,
                     min = 0,
                     method = "simBuyseTest")
        validNumeric(shape.Censoring.C,
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
        lava::distribution(modelT, name0[iterE]) <- lava::weibull.lvm(scale = scale.T[iterE], shape = 1/shape.T[iterE])
        lava::distribution(modelT, nameC[iterE]) <- lava::weibull.lvm(scale = scale.Censoring.T[iterE], shape = 1/shape.Censoring.T[iterE])
        if(CR){
            lava::distribution(modelT, nameCR[iterE]) <- lava::weibull.lvm(scale = scale.CR[iterE], shape = 1/shape.CR[iterE])
            txtSurv <- paste0(name[iterE], "~min(",nameCR[iterE],"=2,",name0[iterE],"=1,",nameC[iterE],"=0)")
        }else{
            txtSurv <- paste0(name[iterE], "~min(",name0[iterE],"=1,",nameC[iterE],"=0)")
        }        
        modelT <- lava::eventTime(modelT, stats::as.formula(txtSurv), nameCensoring[iterE])

        lava::distribution(modelC, name0[iterE]) <- lava::weibull.lvm(scale = scale.C[iterE], shape = 1/shape.C[iterE])
        lava::distribution(modelC, nameC[iterE]) <- lava::weibull.lvm(scale = scale.Censoring.C[iterE], shape = 1/shape.Censoring.C[iterE])
        if(CR){
            lava::distribution(modelC, nameCR[iterE]) <- lava::weibull.lvm(scale = scale.CR[iterE], shape = 1/shape.CR[iterE])
            txtSurv <- paste0(name[iterE], "~min(",nameCR[iterE],"=2,",name0[iterE],"=1,",nameC[iterE],"=0)")
        }else{
            txtSurv <- paste0(name[iterE], "~min(",name0[iterE],"=1,",nameC[iterE],"=0)")
        }
        modelC <- lava::eventTime(modelC, stats::as.formula(txtSurv), nameCensoring[iterE])

        if(CR){
            formula.latent <- as.formula(paste0("~",name0[iterE],"+",nameC[iterE],"+",nameCR[iterE]))
        }else{
            formula.latent <- as.formula(paste0("~",name0[iterE],"+",nameC[iterE]))
        }
        latent(modelT) <- formula.latent
        latent(modelC) <- formula.latent
    }

    ## ** export
    return(list(modelT = modelT, modelC = modelC, latent0 = name0, latentC = nameC, scale.T = scale.T, scale.C = scale.C, shape.T = shape.T, shape.C = shape.C))
}
