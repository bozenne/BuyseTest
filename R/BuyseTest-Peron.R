### BuyseTest-Peron.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 12 2020 (11:10) 
## Version: 
## Last-Updated: nov 18 2020 (14:10) 
##           By: Brice Ozenne
##     Update #: 204
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calcPeron
#' @rdname internal-initialization
calcPeron <- function(data,
                      model.tte, fitter, args,
                      method.score,
                      treatment,
                      level.treatment,
                      endpoint,
                      endpoint.TTE,
                      endpoint.UTTE,
                      status,
                      status.TTE,
                      status.UTTE,
                      D.TTE,
                      D.UTTE,
                      type,
                      strata,
                      threshold,
                      n.strata,
                      precompute,
                      iidNuisance,
                      out){

    ## ** prepare (necessary to compute when resampling)
    if(n.strata == 1){
        ls.indexC <- list(which(data[[treatment]]==0))
        ls.indexT <- list(which(data[[treatment]]==1))
    }else{
        indexC <- which(data[[treatment]]==0)
        indexT <- which(data[[treatment]]==1)
        ls.indexC <- vector(mode = "list", length = n.strata)
        ls.indexT <- vector(mode = "list", length = n.strata)
        for(iStrata in 1:n.strata){
            iIndex.strata <- which(data[["..strata.."]]==iStrata)
            ls.indexC[[iStrata]] <- intersect(indexC,iIndex.strata)
            ls.indexT[[iStrata]] <- intersect(indexT,iIndex.strata)
        }
    }            
    nStrata.control <- lapply(ls.indexC,length)
    nStrata.treatment <- lapply(ls.indexT,length)
    zeroPlus <- 1e-12

    ## ** estimate cumulative incidence function (survival case or competing risk case)
    ls.indexAssociatedEndpoint <- setNames(vector(mode = "list", length = D.UTTE), endpoint.UTTE)
    test.CR <- setNames(vector(mode = "logical", length = D.UTTE), endpoint.UTTE)

    ## prepare formula
    if(is.null(model.tte)){        
        model.tte <- vector(length = D.UTTE, mode = "list")
        names(model.tte) <- endpoint.UTTE
        tofit <- TRUE
        if(length(args)==0){args <- NULL}
        
        if(any(fitter=="prodlim")){
            txt.modelUTTE <- paste0("prodlim::Hist(",endpoint.UTTE,",",status.UTTE,") ~ ",treatment," + ..strata..")
        }else if(any(fitter=="survreg")){
            txt.modelUTTE <- paste0("survival::Surv(",endpoint.UTTE,",",status.UTTE,") ~ ",treatment," + ..strata..")
        }
    }else{
        tofit <- FALSE
    }
    
    ## fit survival model and prepare for extracting survival
    for(iUTTE in 1:D.UTTE){ ## iUTTE <- 1
        ls.indexAssociatedEndpoint[[iUTTE]] <-  which(endpoint == endpoint.UTTE[iUTTE])
        test.CR[iUTTE] <- any(method.score[ls.indexAssociatedEndpoint[[iUTTE]]]==5)

        if(fitter[iUTTE]=="prodlim"){

            if(tofit){
                model.tte[[iUTTE]] <- do.call(prodlim::prodlim, args = c(list(as.formula(txt.modelUTTE[iUTTE]), data = data), args))
            }

        }else if(fitter[iUTTE]=="survreg"){

            if(tofit){
                model.tte[[iUTTE]] <- do.call(survreg::survfit, args = c(list(as.formula(txt.modelUTTE[iUTTE]), data = data), args))
            }
            
        }
        model.tte[[iUTTE]] <- BuyseTTEM(model.tte[[iUTTE]], treatment = treatment, level.treatment = level.treatment, iid = iidNuisance)
    }

    ## ** estimate quantities for scoring pairs
    for(iUTTE in 1:D.UTTE){ ## iUTTE <- 1

        
        ## fitted survival at event timepoints
        for(iS in 1:n.strata){
            for(iT in 0:1){

                iStore <- switch(as.character(iT),
                                 "0" = "survTimeC",
                                 "1" = "survTimeT")

                iTime <- data[list(iT,iS),.SD[[endpoint.UTTE[iUTTE]]],on=c(treatment,"..strata..")]

                iPred.C <- predictWithIndex(model.tte[[iUTTE]], time = iTime, treatment = level.treatment[1], strata = iS)
                iPred.T <- predictWithIndex(model.tte[[iUTTE]], time = iTime, treatment = level.treatment[2], strata = iS)

                
                for(iE in ls.indexAssociatedEndpoint[[iUTTE]]){
                    iPred.C.beforeTau <- predictWithIndex(model.tte[[iUTTE]], time = iTime - threshold[iE], treatment = level.treatment[1], strata = iS)
                    iPred.C.afterTau <- predictWithIndex(model.tte[[iUTTE]], time = iTime - threshold[iE], treatment = level.treatment[1], strata = iS)

                    iPred.T.beforeTau <- predictWithIndex(model.tte[[iUTTE]], time = iTime + threshold[iE], treatment = level.treatment[2], strata = iS)
                    iPred.T.afterTau <- predictWithIndex(model.tte[[iUTTE]], time = iTime + threshold[iE], treatment = level.treatment[2], strata = iS)
                    browser()
                    ## fix pb size + add iid
                    out[[iStore]][[iE]][[iS]] <- cbind("time" = iTime,
                                                       "survivalC-threshold" = iPred.C.beforeTau$survival,
                                                       "survivalC_0" = iPred.C$survival,
                                                       "survivalC+threshold" = iPred.C.afterTau$survival,
                                                       "survivalT-threshold" = iPred.T.beforeTau$survival,
                                                       "survivalT_0" = iPred.T$survival,
                                                       "survivalT+threshold" = iPred.T.afterTau$survival,
                                                       "index.survivalC-threshold" = iPred.C.beforeTau$index,
                                                       "index.survivalC_0" = iPred.C$index,
                                                       "index.survivalC+threshold" = iPred.C.afterTau$index,
                                                       "index.survivalT-threshold" = iPred.T.beforeTau$index,
                                                       "index.survivalT_0" = iPred.T$index,
                                                       "index.survivalT+threshold" = iPred.T.afterTau$index
                                                       )
                }

                
                browser()
                data[list()]
                browser()
            }
        }
        
        
        out[c("survTimeC","survTimeT")] <- calcPeronUTTE(model.tte[[iEndpoint.UTTE]],
                                                         index.start = ls.Mindex.start[[iUTTE]],
                                                         index.stop = ls.Mindex.stop[[iUTTE]],
                                                         test.CR = test.CR[iEndpoint.UTTE],
                                                         lastSurv = out$lastSurv[[ls.indexAssociatedEndpoint[[iUTTE]][1]]],
                                                         level.treatment = level.treatment,
                                                         endpoint = ls.indexAssociatedEndpoint[[iUTTE]],
                                                         threshold = threshold,
                                                         n.strata = n.strata,
                                                         iid = iidNuisance)
        
        ## Treatment group
    }        

    ## ** integrals 
    for(iEndpoint.UTTE in 1:D.UTTE){ ## iEndpoint.TTE <- 1
        ## Control group

        ## Treatment group
    }

    ## ** iid
    
    ## ** export
    return(out)
    
}




######################################################################
### BuyseTest-Peron.R ends here
