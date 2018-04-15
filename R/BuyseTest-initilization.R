## * Documentation initialization functions called by BuyseTest

#' @name internal-intilisation
#' @rdname internal-intilisation
#' @title internal functions for BuyseTest - intilisation
#' @aliases initCensoring initData initN initHypothesis initSpace initStrata initSurvival initThreshold initWscheme
#' 
#' @description Functions called by \code{\link{BuyseTest}} to initialize the arguments.
#' 
#' @details
#' All functions performs test to check the validity of the argument. In addition: \cr
#' 
#' 
#' \code{initCensoring}: If no TTE and censoring is full of NA, set it to NULL. Else when censoring contains NA for the non-TTE endpoints, remove them. \cr
#' \code{initData}: Check that TTE outcome and censoring variable do not contain NAs,
#' display a warning if the other outcomes contain NAs,
#' check that binary/continuous/TTE are indeed what binary/numeric/positive,
#' define design matrix for each group and additional information (e.g. survival). \cr
#' \code{initFormula}: If a formula is given in argument, extract the type of endpoints, the name of the endpoints, the thresholds and the name of the status variables. \cr
#' \code{initStrata}: merge the strata into one with the interaction function and compute the index of the observations relative to each strata. \cr
#' \code{initThreshold}: set default threshold to 1e-12 expect for binary variable where it is set to 1. 
#' Also set null threshold to a small positive value (0->1e-12) 
#' - the rational being we consider a pair favorable if X>Y ie X>=Y+1e-12. When using a threshold e.g. 5 we want X>=Y+5 and not X>Y+5, especially when the measurement is discrete. \cr
#' 
#' 
#' 
#' @keywords function internal BuyseTest

## Functions called by BuyseTest for the initialization

## * applyOperator
#' @rdname internal-intilisation
applyOperator <- function(data, operator, type, endpoint, D){

    validCharacter(operator,
                   valid.values = c("<0",">0"),
                   valid.length = D,
                   method = "BuyseTest")

    n.operatorPerEndpoint <- tapply(operator, endpoint, function(x){length(unique(x))})
    if(any(n.operatorPerEndpoint>1)){
        stop("Cannot have different operator for the same endpoint used at different priorities \n")
    }
    
    operator.endpoint <- setNames(operator, endpoint)[!duplicated(endpoint)]
    name.negative <- names(operator.endpoint)[operator.endpoint=="<0"]
    if(length(name.negative)>0){
        name.negative.binary <- intersect(name.negative, endpoint[type==1])
        if(length(name.negative.binary)>0){
            data[,name.negative.binary] <- -data[,name.negative.binary]+1
        }

        name.negative.other <- setdiff(name.negative, name.negative.binary)
        if(length(name.negative.other)){
            data[,name.negative.other] <- -data[,name.negative.other]
        }
    }
        
    return(data)
}

## * initCensoring
#' @rdname internal-intilisation
initCensoring <- function(censoring,endpoint,type,D,D.TTE,
                          treatment,strata){

  if(D.TTE==0){
    
    if(!is.null(censoring) && all(is.na(censoring))){
      censoring <- NULL
    }
    
    if(!is.null(censoring)){
      stop("BuyseTest : \'censoring\' must be NULL when there are no TTE endpoints \n",
           "propose value : ",paste(censoring,collapse=" "),"\n")
    }
    
  }else{
    
    if(length(censoring) == D){
      
        if(any(!is.na(censoring[type!=3]))){
            stop("BuyseTest : wrong specification of \'censoring\' \n",
                 "\'censoring\' must be NA for binary or continuous endpoints \n",
                 "binary or continuous endoints : ",paste(endpoint[type!=3],collapse=" "),"\n",
                 "proposed \'censoring\' for these endoints : ",paste(censoring[type!=3],collapse=" "),"\n")
        }
      
        if(any(is.na(censoring[type==3])) ){
            stop("BuyseTest : wrong specification of \'censoring\' \n",
                 "\'censoring\' must indicate a variable in data for TTE endpoints \n",
                 "TTE endoints : ",paste(endpoint[type==3],collapse=" "),"\n",
                 "proposed \'censoring\' for these endoints : ",paste(censoring[type==3],collapse=" "),"\n")
        }
      
        censoring <- censoring[type==3] # from now, censoring contains for each time to event endpoint the name of variable indicating censoring (0) or event (1)
    }else{
        
        if(length(censoring) != D.TTE){
            stop("BuyseTest : \'censoring\' does not match \'endpoint\' size \n",
                 "length(censoring) : ",length(censoring),"\n",
                 "length(endpoint) : ",D,"\n")
      
        }

        if(any(is.na(censoring)) ){
            stop("BuyseTest : wrong specification of \'censoring\' \n",
                 "\'censoring\' must indicate a variable in data for TTE endpoints \n",
                 "TTE endoints : ",paste(endpoint[type==3],collapse=" "),"\n",
                 "proposed \'censoring\' for these endoints : ",paste(censoring,collapse=" "),"\n")
        }
        
    }
  }
  
  ## ** export
  return(censoring)
}

## * initData
#' @rdname internal-intilisation
initData <- function(dataT, dataC, type, endpoint, D, censoring, operator,
                     index.strataT, index.strataC, n.strata,          
                     method, D.TTE, threshold, Wscheme = NULL,
                     trace, test = TRUE){
  
  ## ** check NA
  if(test && any(type==3)){
    test.naT <- colSums(is.na(dataT[,c(censoring,endpoint[type==3]),with=FALSE]))
    test.naC <- colSums(is.na(dataC[,c(censoring,endpoint[type==3]),with=FALSE]))
    test.na <- test.naT+test.naC
    if(any(test.na>0)){ 
      stop("BuyseTest : TTE endpoint / censoring variables must not contain NA \n",
           "number of NA in \'data\' : ",sum(test.na>0),"\n",
           "columns with \'data\' : ",paste(c(endpoint,censoring[type==3])[test.na>0],collapse=" "),"\n")
    }
  }
  
    if(test && trace>0 && any(type!=3)){
        test.naT <- colSums(is.na(dataT[,endpoint[type!=3],with=FALSE]))
        test.naC <- colSums(is.na(dataC[,endpoint[type!=3],with=FALSE]))
        test.na <- test.naT+test.naC
        if(any(test.na>0)){
            vec.print <- apply(data.frame(names(test.na),as.double(test.na)),1,paste,collapse=" : ")
            warning("BuyseTest : some binary or continuous endpoints contains NA (variable : number of NA) \n",
                    paste(vec.print, collapse = " \n "))
        }
    
  }

    ## ** endpoint checking : binary type
    indexY <- which(type==1)
    if(test && length(indexY)>0){
        for(iterY in indexY){ ## iterY <- 1
            validNumeric(dataT[[endpoint[iterY]]],
                         name1 = endpoint[iterY],
                         valid.values = 0:1,
                         refuse.NA =  FALSE,
                         valid.length = NULL,
                         method = "BuyseTest")
            validNumeric(dataC[[endpoint[iterY]]],
                         name1 = endpoint[iterY],
                         valid.values = 0:1,
                         refuse.NA =  FALSE,
                         valid.length = NULL,
                         method = "BuyseTest")
        }
    }
  
  ## ** endpoint checking : continuous type
  indexY <- which(type==2)
  if(test && length(indexY)>0){
    for(iterY in indexY){
      validNumeric(dataT[[endpoint[iterY]]],
                   name1 = endpoint[iterY],
                   valid.length = NULL,
                   refuse.NA =  FALSE,
                   method = "BuyseTest")
      validNumeric(dataC[[endpoint[iterY]]],
                   name1 = endpoint[iterY],
                   valid.length = NULL,
                   refuse.NA =  FALSE,
                   method = "BuyseTest")
    }
  }
  
  ## ** endpoint checking : time to event
  indexY <- which(type==3)
  if(test && length(indexY)>0){
    for(iterY in indexY){
      validNumeric(dataT[[endpoint[iterY]]],
                   name1 = endpoint[iterY],
                   valid.length = NULL,
                   method = "BuyseTest")
      validNumeric(dataC[[endpoint[iterY]]],
                   name1 = endpoint[iterY],
                   valid.length = NULL,
                   method = "BuyseTest")
      validNumeric(dataT[[censoring[which(indexY == iterY)]]],
                   name1 = censoring[which(indexY == iterY)],
                   valid.values = c(0,1),
                   valid.length = NULL,
                   method = "BuyseTest")
      validNumeric(dataC[[censoring[which(indexY == iterY)]]],
                   name1 = censoring[which(indexY == iterY)],
                   valid.values = c(0,1),
                   valid.length = NULL,
                   method = "BuyseTest")
    }
  }
  
  ## ** transformation
  M.Treatment <- as.matrix(dataT[,endpoint,with=FALSE]) # matrix of endpoints for the treatment arm 
  M.Control <- as.matrix(dataC[,endpoint,with=FALSE]) # matrix of endpoints for the control arm
  
  if(!is.null(censoring)){
    M.delta_Treatment <- as.matrix(dataT[,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
    M.delta_Control <- as.matrix(dataC[,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
    
    if(method==2){ # "Efron"
      for(iter_strata in 1:n.strata){
        
        Mstrata.Treatment <- M.Treatment[index.strataT[[iter_strata]]+1,,drop=FALSE]
        Mstrata.Control <- M.Control[index.strataC[[iter_strata]]+1,,drop=FALSE]
        
        # set last observation for each TTE endpoint to non-censored
        for(iter_endpointTTE in 1:D.TTE){
          indexT_maxCensored <- which(Mstrata.Treatment[,which(type==3)[iter_endpointTTE]]==max(Mstrata.Treatment[,which(type==3)[iter_endpointTTE]]))
          M.delta_Treatment[index.strataT[[iter_strata]][indexT_maxCensored]+1,iter_endpointTTE] <- 1
          indexC_maxCensored <- which(Mstrata.Control[,which(type==3)[iter_endpointTTE]]==max(Mstrata.Control[,which(type==3)[iter_endpointTTE]]))
          M.delta_Control[index.strataC[[iter_strata]][indexC_maxCensored]+1,iter_endpointTTE] <- 1
        }
      }
    }
    
  }else{ # if the is no time to event variables
    M.delta_Treatment <- matrix(-1,1,1) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
    M.delta_Control <- matrix(-1,1,1) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
  }
  
  ## ** KM imputation
  if(method %in% 1:3){# c("Peto","Efron","Peron")
    
    endpoint.TTE <- endpoint[type==3] # vector of variable names of the TTE endpoints
    
    ## *** design matrix for the weights
    if(is.null(Wscheme)){
        res_init <- initWscheme(D=D,
                                endpoint=endpoint,
                                endpoint.TTE=endpoint.TTE,
                                D.TTE=D.TTE,
                                threshold=threshold,
                                type=type)
      Wscheme <- res_init$Wscheme  
      index_survivalM1 <- res_init$index_survivalM1
      threshold_TTEM1 <- res_init$threshold_TTEM1       
    }
    ## *** Survival estimate using Kaplan Meier    
    res_init <- initSurvival(M.Treatment=M.Treatment,
                             M.Control=M.Control,
                             M.delta_Treatment=M.delta_Treatment,
                             M.delta_Control=M.delta_Control,
                             endpoint=endpoint,
                             D.TTE=D.TTE,
                             type=type,
                             threshold=threshold,
                             index.strataT=index.strataT,
                             index.strataC=index.strataC,
                             n.strata=n.strata,   
                             method=method)
    
    list_survivalT <- res_init$list_survivalT
    list_survivalC <- res_init$list_survivalC
    
  }else{
    Wscheme <- matrix()  # factice design matrix for the weights. Will be sent to the C++ arguments to fill the argument but not used by the function.
    list_survivalT <- list() # factice list. Will be sent to the C++ arguments to fill the argument but not used by the function.
    list_survivalC <- list() # factice list. Will be sent to the C++ arguments to fill the argument but not used by the function.
    index_survivalM1 <- numeric(0) # factice vector. Will be sent to the C++ arguments to fill the argument but not used by the function.
    threshold_TTEM1 <- numeric(0)  # factice vector. Will be sent to the C++ arguments to fill the argument but not used by the function.
    }
  
  ## ** export 
  res <- list()
  res$M.Treatment <-  M.Treatment
  res$M.Control <-  M.Control
  res$M.delta_Treatment <-  M.delta_Treatment
  res$M.delta_Control <-  M.delta_Control
  res$Wscheme <- Wscheme
  res$threshold_TTEM1 <- threshold_TTEM1
  res$index_survivalM1 <- index_survivalM1
  res$list_survivalT <- list_survivalT
  res$list_survivalC <- list_survivalC
  
  return(res)
}

## * initWscheme
#' @rdname internal-intilisation
initWscheme <- function(endpoint,D,endpoint.TTE,D.TTE,threshold,type){

  if(D>1){
        
    Wscheme <- matrix(NA,nrow=D.TTE,ncol=D-1) # design matrix indicating to combine the weights obtained at differents TTE endpoints
    rownames(Wscheme) <- paste("weigth of ",endpoint.TTE,"(",threshold[type==3],")",sep="")
    colnames(Wscheme) <- paste("for ",endpoint[-1],"(",threshold[-1],")",sep="")
    
    index_survivalM1 <- rep(-1,D.TTE) # index of previous TTE endpoint (-1 if no previous TTE endpoint i.e. endpoint has not already been used)
    threshold_TTEM1 <- rep(-1,D.TTE) # previous threshold (-1 if no previous threshold i.e. endpoint has not already been used)
    
    iter_endpoint.TTE <- if(type[1]==3){1}else{0} #  index_survivalM1 and  threshold_TTEM1 are -1 even if the first endoint is a survival endpoint
    
    for(iter_endpoint in 2:D){   
      
      if(type[iter_endpoint]==3){ 
        iter_endpoint.TTE <- iter_endpoint.TTE + 1
      }
      
      # select valid rows
      index_rowTTE <- which(paste(endpoint.TTE,threshold[type==3],sep="_") %in% paste(endpoint,threshold,sep="_")[1:(iter_endpoint-1)])
      if(length(index_rowTTE)==0){next} # not yet TTE endpoints (no valid rows)
      Wscheme[index_rowTTE,iter_endpoint-1] <- 0 # potential weights
      
      # keep only the last repeated endpoints
      index_rowTTE <- sapply(unique(endpoint.TTE[index_rowTTE]),
                                   function(x){utils::tail(index_rowTTE[endpoint.TTE[index_rowTTE]==x],1)})
      
      # if survival endpoint remove similar endpoints
      if(endpoint[iter_endpoint] %in% endpoint.TTE[index_rowTTE]){
        index_survivalM1[iter_endpoint.TTE] <- index_rowTTE[endpoint.TTE[index_rowTTE]==endpoint[iter_endpoint]]-1
        threshold_TTEM1[iter_endpoint.TTE] <- threshold[type==3][index_survivalM1[iter_endpoint.TTE]+1]
        index_rowTTE <- setdiff(index_rowTTE,index_survivalM1[iter_endpoint.TTE]+1)
      }
      
      # update Wscheme
      if(length(index_rowTTE)>0){
        Wscheme[index_rowTTE,iter_endpoint-1] <- 1
      }
      
    }

  }else{
    Wscheme <- matrix(nrow=0,ncol=0)
    index_survivalM1 <- numeric(0)
    threshold_TTEM1 <- numeric(0)
  }
  
  ## ** export
  res <- list()
  res$Wscheme <- Wscheme
  res$index_survivalM1 <- index_survivalM1
  res$threshold_TTEM1 <- threshold_TTEM1
  
  return(res)
  
}

## * initFormula
#' @rdname internal-intilisation
#' @examples 
#' BuyseTest:::initFormula(Treatment~B(var1)+C(var2,10)+T(var3,15,exa))
#' BuyseTest:::initFormula(Treatment~ grade + Bin(var1)+C(var2,10)+T(var3,15,exa))
#' BuyseTest:::initFormula(Treatment~ grade + Bin(var1)+C(var2,10)+T(var3,15,exa) + treatment)
initFormula <- function(x){

  validClass(x, valid.class = "formula")
    
  ## ** extract treatment
  treatment <- setdiff(all.vars(x), all.vars(stats::delete.response(stats::terms(x))))
  if(length(treatment)!=1){
    stop("initFormula: there must be exactly one response variable in formula\n",
         "number of response variables founded: ",length(treatment),"\n")
  }
  
  if(length(as.character(x))!=3){
    stop("initFormula: formula with unexpected length, as.character(x) should have length 3\n",
         "length founded: ",length(as.character(x)),"\n")
  }
  
    ## ** restrict to the right side of the formula
    x.rhs <- as.character(x)[3]
  
    ## remove all blancks
    x.rhs <- gsub("[[:blank:]]", "", x.rhs)
    vec.x.rhs <- unlist(strsplit(x.rhs, split = "+", fixed = TRUE))

    ## find all element in the vector corresponding to endpoints (i.e. ...(...) )
    ## \\w* any letter/number
    ## [[:print:]]* any letter/number/punctuation/space
    index.endpoint <- grep("^\\w*\\([[:print:]]*\\)$", vec.x.rhs)
    index.strata <- setdiff(1:length(vec.x.rhs), index.endpoint)

    ## ** strata variables
    if(length(index.strata)==0){
        strata <- NULL
    }else{
        strata <- vec.x.rhs[index.strata]
    }
    
    ## ** number of endpoint variables    
    vec.x.endpoint <- vec.x.rhs[index.endpoint]
    n.endpoint <- length(vec.x.endpoint)
    if(length(n.endpoint)==0){
        stop("initFormula: x must contains endpoints \n",
             "nothing of the form type(endpoint,threshold,censoring) found in the formula \n")
    }

    ## ** extract endpoints and additional arguments 
    threshold <- rep(NA, n.endpoint)
    censoring <- rep(NA, n.endpoint)
    endpoint <- rep(NA, n.endpoint)
    operator <- rep(">0", n.endpoint)
    validArgs <- c("endpoint","threshold","censoring","operator")

    ## split around parentheses
    ls.x.endpoint <- strsplit(vec.x.endpoint, split = "(", fixed = TRUE)

    type <- character(length = n.endpoint)
    for(iE in 1:n.endpoint){
        ## extract type
        type[iE] <- ls.x.endpoint[[iE]][1]
        
        ## get each argument
        iVec.args <- strsplit(gsub(")", replacement = "",ls.x.endpoint[[iE]][2]),
                              split = ",", fixed = TRUE)[[1]]
        n.args <- length(iVec.args)
        
        ## check size
        if(n.args==0){
            stop("initFormula: invalid formula \n",
                 vec.x.rhs[iE]," must contain the name of the endpoint between the parentheses \n"
                 )
        }
        if(n.args>3){
            stop("initFormula: invalid formula \n",
                 x[iE]," has too many arguments (maximum 4: endpoint, threshold, censoring variable, operator) \n")
        }

        ## extract name of each argument
        iIndex.name <- grep("=",iVec.args)
        iArg <- gsub("^[[:print:]]*=", replacement = "", iVec.args)        
        iName <- rep(as.character(NA),n.args)
        
        ## use existing names
        if(length(iIndex.name)>0){
            iiName <- gsub("=[[:print:]]*$","",iVec.args[iIndex.name])
            iName[iIndex.name] <- iiName
            if(any(iiName %in% validArgs == FALSE)){
                stop("initFormula: invalid formula \n",
                     vec.x.rhs[iE]," contains arguments that are not endpoint, threshold, censoring \n")
            }
            if( any(duplicated(iiName)) ){
                stop("initFormula: invalid formula \n",
                     vec.x.rhs[iE]," contains arguments with the same name \n")
            }
        }else{
            iiName <- NULL
        }
        
        ## add missing names
        n.missingNames <- n.args - length(iiName) 
        if(n.missingNames>0){
            iName[setdiff(1:n.args,iIndex.name)] <- setdiff(validArgs,iiName)[1:n.missingNames]
        }

        ## extract arguments
        endpoint[iE] <- gsub("\"","",iArg[iName=="endpoint"])
        if("threshold" %in% iName){
            threshold[iE] <- as.numeric(iArg[iName=="threshold"])
        }
        if("censoring" %in% iName){
            censoring[iE] <- gsub("\"","",iArg[iName=="censoring"])
        }
        if("operator" %in% iName){
            operator[iE] <- gsub("\"","",iArg[iName=="operator"])
        }
    }

    ## ** export
    return(list(treatment = treatment,
                type = type,
                endpoint = endpoint,
                threshold = threshold,
                censoring = censoring,
                operator = operator,
                strata = strata))
}

## * initStrata
#' @rdname internal-intilisation
initStrata <- function(strata,
                       dataT,dataC,n.Treatment,n.Control,
                       endpoint,censoring){
  
  if(!is.null(strata)){
    strataT <- interaction(dataT[,strata,with=FALSE],drop = TRUE,lex.order=FALSE,sep=".") # strata variable for the treatment arm transformed into factor
    strataC <- interaction(dataC[,strata,with=FALSE],drop = TRUE,lex.order=FALSE,sep=".") # strata variable for the control arm transformed into factor
    levels.strata <- levels(strataT) # extraction of the levels of the strata variables
    
    if(length(levels.strata) != length(levels(strataC)) || any(levels.strata != levels(strataC))){
      stop("BuyseTest : wrong specification of \'strata\' \n",
           "different levels between Control and Treatment \n",
           "levels(strataT) : ",paste(levels(strataT),collapse=" "),"\n",
           "levels(strataC) : ",paste(levels(strataC),collapse=" "),"\n")
    }
    
    strataT <- as.numeric(strataT) # strata variable for the treatment arm converted into numeric 
    strataC <- as.numeric(strataC) # strata variable for the control arm converted into numeric
    
  }else{ # if there is no strata variable the same strata is used for all patient
    strataT <- rep(1,n.Treatment) # all patient of the treatment arm are in the same strata : strata 1
    strataC <- rep(1,n.Control)  # all patient of the control arm are in the same strata : strata 1
    levels.strata <- 1 # the only strata is strata 1
  }
  
  n.strata <- length(levels.strata) # number of strata
  index.strataT <- lapply(1:n.strata,function(x){which(x==strataT)-1}) # for each strata, the index of the patients belonging to this strata in the treatment arm is stored in an element of the list. 
  index.strataC <- lapply(1:n.strata,function(x){which(x==strataC)-1}) # for each strata, the index of the patients belonging to this strata in the control arm is stored in an element of the list. 
  # Because of the minus one after the which, both indexes begin at 0. This is compulsory for C++ compatibility
  
  # check there is at least one element per arm in each strata
  if(any(unlist(lapply(index.strataT,length))==0)){
    stop("BuyseTest : some strata contain no patient in the treatment arm \n",
         "strata : ",paste(levels.strata[which(unlist(lapply(index.strataT,length))==0)],collapse=" ")," \n")
  }

  if(any(unlist(lapply(index.strataC,length))==0)){
    stop("BuyseTest : some strata contain no patient in the control arm \n",
         "strata : ",paste(levels.strata[which(unlist(lapply(index.strataC,length))==0)],collapse=" ")," \n")
  }
  
  ## ** export
  res <- list()
  res$index.strataT <- index.strataT
  res$index.strataC <- index.strataC
  res$n.strata <- n.strata
  res$levels.strata <- levels.strata
  
  return(res)
}

## * initThreshold
#' @rdname internal-intilisation
initThreshold <- function(threshold,type,D,endpoint){
  
    ## ** initialize threshold
    if(is.null(threshold)){
        threshold <- rep(10^{-12},D)  # if no treshold is proposed all threshold are by default set to 10^{-12}
        if(any(type==1)){threshold[type==1] <- 1/2} # except for threshold corresponding to binary endpoints that are set to NA.
    }else{
        if(any(is.na(threshold[type==1]))){
            index.tempo <- intersect(which(is.na(threshold)),which(type==1))
            threshold[index.tempo] <- 1/2
        }
        if(any(abs(stats::na.omit(threshold))<10^{-12})){
            threshold[which(abs(threshold)<10^{-12})] <- 10^{-12}
        }
        if(any(is.na(threshold))){
            threshold[is.na(threshold)] <- 10^{-12}
        }
        validNumeric(threshold,
                     valid.length = D,
                     min = 0,
                     refuse.NA = FALSE ,
                     method = "BuyseTest")
    }
    
    ## ** Only accept hreshold 1/2 for binary outcomes
    if(any(threshold[type==1]!=1/2)){
        stop("BuyseTest : wrong specification of \'threshold\' \n",
             "\'threshold\' must be NA for binary endpoints (or equivalently 1/2) \n",
             "proposed \'threshold\' : ",paste(threshold[type==1],collapse=" "),"\n",
             "binary endpoint(s) : ",paste(endpoint[type==1],collapse=" "),"\n")
    }
  
    ## ** Check that the thresholds related to the same endoints are strictly decreasing
    ## is.unsorted(rev(2:1))
    ## is.unsorted(rev(1:2))

    vec.test <- tapply(threshold,endpoint, function(x){
        test.unsorted <- is.unsorted(rev(x))
        test.duplicated <- any(duplicated(x))
        return(test.unsorted+test.duplicated)
    })
            
    if(any(vec.test>0)){   
        stop("BuyseTest : wrong specification of \'endpoint\' or \'threshold\' \n",
             "endpoints must be used with strictly decreasing threshold when re-used with lower priority \n",
             "problematic endpoints: \"",paste0(names(vec.test)[vec.test>0], collapse = "\" \""),"\"\n")        
    }

    
  ## ** export
  return(threshold)
}

## * initSurvival
#' @rdname internal-intilisation
#' @export
initSurvival <- function(M.Treatment,M.Control,M.delta_Treatment,M.delta_Control,
                         endpoint,D.TTE,type,threshold,
                         index.strataT,index.strataC,n.strata, 
                         method){
  
  ## ** conversion to R index
  index.strataT <- lapply(index.strataT,function(x){x+1})
  index.strataC <- lapply(index.strataC,function(x){x+1})
  
  ## ** preparing 
  if(method %in% 2:3){ # c("Efron","Peron")
    colnamesT <- c(paste("SurvivalT_TimeT",c("-threshold","_0","+threshold"),sep=""),
                   paste("SurvivalC_TimeT",c("-threshold","_0","+threshold"),sep=""),
                   "SurvivalT_TimeT_0(ordered)","SurvivalC_TimeT-threshold(ordered)","SurvivalC_TimeT+threshold(ordered)",
                   "time_control(ordered)","events(ordered)") 
    
    colnamesC <- c(paste("SurvivalC_TimeC",c("-threshold","_0","+threshold"),sep=""),
                   paste("SurvivalT_TimeC",c("-threshold","_0","+threshold"),sep=""),
                   "SurvivalC_TimeC_0(ordered)","SurvivalT_TimeC-threshold(ordered)","SurvivalT_TimeC+threshold(ordered)",
                   "time_control(ordered)","events(ordered)") 
    
    list_survivalT <- rep(list(matrix(NA,nrow=nrow(M.Treatment),ncol=11,dimnames=list(1:nrow(M.Treatment),colnamesT))),times=D.TTE) # list of matrix of survival data for each endpoint in the treatment arm
    list_survivalC <- rep(list(matrix(NA,nrow=nrow(M.Control),ncol=11,dimnames=list(1:nrow(M.Control),colnamesC))),times=D.TTE) # list of matrix of survival data for each endpoint in the control arm
  }else if(method == 1){ # method == "Peto"
    colnamesT <- paste("Survival_TimeT",c("-threshold","_0","+threshold"),sep="")
    colnamesC <- paste("Survival_TimeC",c("-threshold","_0","+threshold"),sep="")
    
    list_survivalT <- rep(list(matrix(NA,nrow=nrow(M.Treatment),ncol=3,dimnames=list(1:nrow(M.Treatment),colnamesT))),times=D.TTE) # list of matrix of survival data for each endpoint in the treatment arm
    list_survivalC <- rep(list(matrix(NA,nrow=nrow(M.Control),ncol=3,dimnames=list(1:nrow(M.Control),colnamesC))),times=D.TTE) # list of matrix of survival data for each endpoint in the control arm
  }
  
  ## ** loop 
  for(iter_strata in 1:n.strata){

    Mstrata.Treatment <- M.Treatment[index.strataT[[iter_strata]],,drop=FALSE]
    Mstrata.Control <- M.Control[index.strataC[[iter_strata]],,drop=FALSE]
    Mstrata.delta_Treatment <- M.delta_Treatment[index.strataT[[iter_strata]],,drop=FALSE]
    Mstrata.delta_Control <- M.delta_Control[index.strataC[[iter_strata]],,drop=FALSE]
    
    for(iter_endpointTTE in 1:D.TTE){
     
      if(method %in% 2:3){
   
        ## **** treatment
        df.treatment <- data.frame(time = Mstrata.Treatment[,endpoint[type==3][iter_endpointTTE]], # store the event times of the treatment arm for a given TTE endpoint
                                   status = Mstrata.delta_Treatment[,iter_endpointTTE])
        
        resKMT_tempo <- prodlim::prodlim(prodlim::Hist(time,status)~1, data = df.treatment)
        # prefer sort(time_treatment) to resKMT_tempo$time to avoid rounding error
        
        if(all(Mstrata.delta_Treatment[which(df.treatment$time==max(df.treatment$time)),iter_endpointTTE]==1)){ # step interpolation of the survival function (last = event)
          
          survestKMT_tempo <- stats::approxfun(x=unique(sort(df.treatment$time)),y=resKMT_tempo$surv,
                                               yleft=1,yright=0,f=0, method = "constant") # 0 at infinity if last event is a death      
          
        }else{ # step interpolation of the survival function (last = censoring)
          
          survestKMT_tempo <- stats::approxfun(x=unique(sort(df.treatment$time)),y=resKMT_tempo$surv,
                                               yleft=1,yright=NA,f=0, method = "constant") # NA at infinity if last event is censored
          
        }
        
        ## **** control
        df.control <- data.frame(time = Mstrata.Control[,endpoint[type==3][iter_endpointTTE]], # store the event times of the treatment arm for a given TTE endpoint
                                 status = Mstrata.delta_Control[,iter_endpointTTE])
        
        resKMC_tempo <- prodlim::prodlim(prodlim::Hist(time,status)~1, data = df.control) 
        
        if(all(Mstrata.delta_Control[which(df.control$time==max(df.control$time)),iter_endpointTTE]==1)){  # step interpolation of the survival function (last = event)
          
          survestKMC_tempo <- stats::approxfun(x=unique(sort(df.control$time)), y=resKMC_tempo$surv,
                                               yleft=1,yright=0,f=0, method= "constant") # 0 at infinity if last event is a death
          
        }else{ # step interpolation of the survival function (last = censoring)
          
          survestKMC_tempo <- stats::approxfun(x=unique(sort(df.control$time)),y=resKMC_tempo$surv,
                                               yleft=1,yright=NA,f=0, method= "constant") # NA at infinity if last event is censored      
          
        }
        
        ## *** survival
        ## **** treatment
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],1] <- survestKMT_tempo(df.treatment$time-threshold[type==3][iter_endpointTTE]) # survival at t - tau, tau being the threshold
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],2] <- survestKMT_tempo(df.treatment$time) # survival at t
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],3] <- survestKMT_tempo(df.treatment$time+threshold[type==3][iter_endpointTTE]) # survival at t + tau
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],4] <- survestKMC_tempo(df.treatment$time-threshold[type==3][iter_endpointTTE]) # survival at t - tau
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],5] <- survestKMC_tempo(df.treatment$time) # survival at t
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],6] <- survestKMC_tempo(df.treatment$time+threshold[type==3][iter_endpointTTE]) # survival at t + tau
        
        order_time_treatment <- order(df.treatment$time)
        
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],7:9] <- list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],c(2,4,6),drop=FALSE][order_time_treatment,,drop=FALSE] # St[t_treatment], Sc[t_treatment-threshold,t_treatment+threshold]
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],10] <-  df.treatment$time[order_time_treatment]
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],11] <-   Mstrata.delta_Treatment[order_time_treatment,iter_endpointTTE] # to compute the integral (Efron)
        
        rownames(list_survivalT[[iter_endpointTTE]])[index.strataT[[iter_strata]]] <- df.treatment$time
        
        ## **** control
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],1] <- survestKMC_tempo(df.control$time-threshold[type==3][iter_endpointTTE]) # survival at t - tau
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],2] <- survestKMC_tempo(df.control$time) # survival at t
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],3] <- survestKMC_tempo(df.control$time+threshold[type==3][iter_endpointTTE]) # survival at t + tau
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],4] <- survestKMT_tempo(df.control$time-threshold[type==3][iter_endpointTTE]) # survival at t - tau, tau being the threshold
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],5] <- survestKMT_tempo(df.control$time) # survival at t
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],6] <- survestKMT_tempo(df.control$time+threshold[type==3][iter_endpointTTE]) # survival at t + tau                                              
        
        order_time_control <- order(df.control$time)
        
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],7:9] <- list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],c(2,4,6),drop=FALSE][order_time_control,,drop=FALSE] # Sc[t_control], St[t_control-threshold,t_control+threshold]
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],10] <- df.control$time[order_time_control]
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],11] <- Mstrata.delta_Control[order_time_control,iter_endpointTTE] # to compute the integral (Efron)
        
        rownames(list_survivalC[[iter_endpointTTE]])[index.strataC[[iter_strata]]] <- df.control$time    
    
      }else if(method == 1){ # Peto
        time_treatment <- Mstrata.Treatment[,endpoint[type==3][iter_endpointTTE]] # store the event times of the treatment arm for a given TTE endpoint
        time_control <- Mstrata.Control[,endpoint[type==3][iter_endpointTTE]] # store the event times of the control arm for a given TTE endpoint
        
        # prefer c(time_treatment,time_control) to time_all to avoid rounding error
        df.all <- data.frame(time = c(time_treatment,time_control),
                             status = c(Mstrata.delta_Treatment[,iter_endpointTTE],Mstrata.delta_Control[,iter_endpointTTE]))
        
        resKM_tempo <- prodlim::prodlim(prodlim::Hist(time,status)~1, data = df.all)
        
        if(all(df.all$status[which(df.all$time==max(df.all$time))]==1)){ # step interpolation of the survival function
          survestKM_tempo <- stats::approxfun(x=unique(sort(c(time_treatment,time_control))),y=resKM_tempo$surv,
                                              yleft=1,yright=0,f=0, method= "constant") # 0 at infinity if last event is a death
          
        }else{
          survestKM_tempo <- stats::approxfun(x=unique(sort(c(time_treatment,time_control))),y=resKM_tempo$surv,
                                              yleft=1,yright=NA,f=0, method= "constant") # NA at infinity if last event is a death
        }
        
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],1] <- survestKM_tempo(time_treatment-threshold[type==3][iter_endpointTTE]) # survival at t - tau, tau being the threshold
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],2] <- survestKM_tempo(time_treatment) # survival at t
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],3] <- survestKM_tempo(time_treatment+threshold[type==3][iter_endpointTTE]) # survival at t + tau
       
        rownames(list_survivalT[[iter_endpointTTE]])[index.strataT[[iter_strata]]] <- time_treatment
        
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],1] <- survestKM_tempo(time_control-threshold[type==3][iter_endpointTTE]) # survival at t - tau
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],2] <- survestKM_tempo(time_control) # survival at t
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],3] <- survestKM_tempo(time_control+threshold[type==3][iter_endpointTTE]) # survival at t + tau
        
        rownames(list_survivalC[[iter_endpointTTE]])[index.strataC[[iter_strata]]] <- time_control    
        
      }
      
    }
  }
  
  return(list(list_survivalT=list_survivalT,
              list_survivalC=list_survivalC))
  
}




## * NOT USED
## ** Function initN
#' @rdname internal-intilisation
initN <- function(n){
  
  ## ** test
  
  if(length(n)!=3){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "\'n\' must have length 3 \n",
         "length(n) : ",length(n),"\n")
  }
  
  if(any(n %% 1>0) || any(n<=0)){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "\'n\' must contains strictly positive integers \n",
         "proposed n : ",paste(n,collpase=" "),"\n")
  }
  
  if(is.null(names(n))){names(n) <- c("from","to","by")}
  
  if(any(names(n) %in% c("from","to","by","length.out")==FALSE)){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "elements of \'n\' must be named with \"from\", \"to\", \"by\" or \"length.out\" \n",
         "invalid names : ",paste(names(n)[names(n) %in% c("from","to","by","length.out")==FALSE] ,collpase=" "),"\n")    
  }
  
  if(length(names(n)) != length(unique(names(n)))){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "elements in \'n\' have the same names \n")    
  }
  
  if("from" %in% names(n) ==FALSE){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "\'n\' must contain an element named \"from\" \n",
         "names(n) : ",paste(names(n) ,collpase=" "),"\n")    
  }
  
  ## ** computation
  
  test_to <- "to" %in% names(n)
  test_by <- "by" %in% names(n)
  test_length.out <- "length.out" %in% names(n)
  if( test_to && test_by){
    n <- seq(from=n["from"],to=n["to"],by=n["by"])
  }
  if( test_to && test_length.out){
    n <- seq(from=n["from"],to=n["to"],length.out=n["length.out"])
  }
  if( test_by && test_length.out){
    n <- seq(from=n["from"],by=n["by"],length.out=n["length.out"])
  }
  
  ## ** export
  return(n)
}

## ** Function initHypothesis
#' @rdname internal-intilisation
initHypothesis <- function(hypothesis,type,D){
  
  if(!is.list(hypothesis) && D==1){
    hypothesis <- list(hypothesis)
  }
  
  if(!is.list(hypothesis) || length(hypothesis)!=D){
    stop("BuyseTest : proposed \'hypothesis\' does not match \'type\' \n",
         "\'hypothesis\' must be a list with ",D," elements \n",
         "is(hypothesis) : ",paste(is(hypothesis),collapse=" "),"\n",
         "length(hypothesis) : ",length(hypothesis),"\n"
    )    
  }
  
  valid_parameters <- list(binary=c("proba_t","proba_c"),
                           continuous=c("mu_t","mu_c","sigma"),
                           TTE=c("lambda_t","lambda_c","T_inclusion","T_followUp")
  )
  
  valid_range <- list(proba_t=c(0,1),
                      proba_c=c(0,1),
                      sigma=c(0,Inf),
                      lambda_t=c(0,Inf),
                      lambda_c=c(0,Inf),
                      T_inclusion=c(0,Inf),
                      T_followUp=c(0,Inf))
  
  for(iter_outcome in 1:D){
    
    parameters_tempo <- valid_parameters[[type[iter_outcome]]]
    
    if(length(hypothesis[[iter_outcome]]) != length(parameters_tempo)){
      stop("BuyseTest : proposed \'hypothesis\' does not match \'type\' \n",
           "type of outcome ",iter_outcome," : ",names(valid_parameters)[type[iter_outcome]]," \n",
           "number of parameters requested for this type of outcome : ",length(parameters_tempo)," \n",
           "length(hypothesis[[",iter_outcome,"]]) : ",length(hypothesis[[iter_outcome]])," \n")   
    }
    
    if(is.null(names(hypothesis[[iter_outcome]]))){
      names(hypothesis[[iter_outcome]]) <- parameters_tempo
    }
    
    if(any(parameters_tempo %in% names(hypothesis[[iter_outcome]]) == FALSE)){
      stop("BuyseTest : wrong specification of \'hypothesis\' \n",
           "type of outcome ",iter_outcome," : ",names(valid_parameters)[type[iter_outcome]]," \n",
           "parameters requested for this type of outcome : \"",paste(parameters_tempo,collapse="\" \""),"\" \n",
           "names(hypothesis[[",iter_outcome,"]]) : \"",paste(names(hypothesis[[iter_outcome]])[names(hypothesis[[iter_outcome]]) %in% parameters_tempo == FALSE],collapse="\" \""),"\" \n")   
    }
    
    for(iter_param in 1:length(parameters_tempo)){
      if(parameters_tempo[iter_param] %in% names(valid_range)){
        test.inf <- valid_range[[parameters_tempo[iter_param]]][1]>=hypothesis[[iter_outcome]][parameters_tempo[iter_param]]
        test.sup <- valid_range[[parameters_tempo[iter_param]]][2]<=hypothesis[[iter_outcome]][parameters_tempo[iter_param]]
        if(test.inf || test.sup){
          stop("BuyseTest : unvalid parameter value in \'hypothesis\' \n",
               "possible range of values : ",paste(valid_range[[parameters_tempo[iter_param]]],collapse=" ")," \n",
               "proposed value of hypothesis[[",iter_outcome,"]][",parameters_tempo[iter_param],"] :  ",hypothesis[[iter_outcome]][parameters_tempo[iter_param]]," \n"
          )  
        }
      }
    }
    
    
  }
  
  ### * export
  return(hypothesis)
  
}


## ** Function initSpace
#' @rdname internal-intilisation
initSpace  <- function(nchar){  
  return(sapply(nchar,function(x){if(x>0){do.call(paste,c(as.list(rep(" ",x)),sep=""))}else{""}}))  
}

