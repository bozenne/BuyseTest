### BuyseTest-Peron.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 12 2020 (11:10) 
## Version: 
## Last-Updated: nov 30 2020 (16:48) 
##           By: Brice Ozenne
##     Update #: 288
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

    zeroPlus <- 1e-12

    ## ** fit model for the cumulative incidence function (survival case or competing risk case)
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
        for(iStrata in 1:n.strata){

            for(iTreat in level.treatment){
                iStoreJump <- c("survJumpC","survJumpT")[(iTreat==level.treatment[2])+1]
                iStoreTime <- c("survTimeC","survTimeT")[(iTreat==level.treatment[2])+1]
                iStoreP <- c("p.C","p.T")[(iTreat==level.treatment[2])+1]

                iTreat.num <- as.numeric(iTreat==level.treatment[2])                
                iTime <- data[list(iTreat.num,iStrata),.SD[[endpoint.UTTE[iUTTE]]],on=c(treatment,"..strata..")]
                iPred.C <- predict(model.tte[[iUTTE]], time = iTime, treatment = level.treatment[1], strata = iStrata)
                iPred.T <- predict(model.tte[[iUTTE]], time = iTime, treatment = level.treatment[2], strata = iStrata)
                
                for(iEndpoint in ls.indexAssociatedEndpoint[[iUTTE]]){
                    if(test.CR){
                        
                    }else{
                        ## *** survival at jump time
                        iTime.jump <- model.tte[[iUTTE]]$peron$jumpSurvHaz[[iStrata]][[iTreat]]$time.jump
                        
                        iSurvTau.jump <- predict(model.tte[[iUTTE]], time = iTime.jump + threshold[iEndpoint], treatment = setdiff(level.treatment, iTreat),
                                                 strata = iStrata, iid = iidNuisance)
                        iSurvBefore.jump <- predict(model.tte[[iUTTE]], time = iTime.jump-zeroPlus, treatment = iTreat,
                                                    strata = iStrata, iid = iidNuisance)
                        iSurvAfter.jump <- predict(model.tte[[iUTTE]], time = iTime.jump+zeroPlus, treatment = iTreat,
                                                   strata = iStrata, iid = iidNuisance)

                        if(iidNuisance){
                            out$iid[[iStoreJump]][[iUTTE]][[iStrata]] <- model.tte[[iUTTE]]$peron$survival.iid[[iStrata]][[iTreat]]
                            out[[iStoreP]][iStrata, iEndpoint] <- NCOL(out$iid[[iStoreJump]][[iUTTE]][[iStrata]])
                        }

                        ## rm after last event if NA and keep last survival value
                        out$lastSurv[[iEndpoint]][iStrata,iTreat.num+1] <- 1-model.tte[[iUTTE]]$peron$last.estimate[iStrata,which(level.treatment==iTreat)] ## 1- because it is cif not survival that is output

                        ## export
                        out[[iStoreJump]][[iEndpoint]][[iStrata]] <- cbind(time = iTime.jump, ## jump time
                                                                           survival = iSurvTau.jump$survival, 
                                                                           dSurvival = iSurvAfter.jump$survival - iSurvBefore.jump$survival,
                                                                           index.survival = iSurvTau.jump$index, ## index of the survival parameter at t+\tau
                                                                           index.dsurvival1 = iSurvBefore.jump$index, ## index of the survival parameter before the jump
                                                                           index.dsurvival2 = iSurvAfter.jump$index) ## index of the survival parameter after the jump

                        ## *** survival at observation time (+/- threshold)
                        iPred.C.beforeTau <- predict(model.tte[[iUTTE]], time = iTime - threshold[iEndpoint], treatment = level.treatment[1], strata = iStrata)
                        iPred.C.afterTau <- predict(model.tte[[iUTTE]], time = iTime + threshold[iEndpoint], treatment = level.treatment[1], strata = iStrata)

                        iPred.T.beforeTau <- predict(model.tte[[iUTTE]], time = iTime - threshold[iEndpoint], treatment = level.treatment[2], strata = iStrata)
                        iPred.T.afterTau <- predict(model.tte[[iUTTE]], time = iTime + threshold[iEndpoint], treatment = level.treatment[2], strata = iStrata)

                        out[[iStoreTime]][[iEndpoint]][[iStrata]] <- cbind("time" = iTime,
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
                    } ## End if CR
                } ## End if iEndpoint
            } ## End if treatment
        } ## End if strata
    } ## End if UTTE
    
    ## ** integrals 
    for(iEndpoint in 1:length(endpoint)){ ## iEndpoint <- 1

        if(!precompute || method.score[iEndpoint]!=4){next} ## only for survival - not (yet!) available for the competing risk case
    
        for(iStrata in 1:n.strata){  ## iStrata <- 1

            ## compute integral at any jump time 
            ls.intC <- calcIntegralSurv2_cpp(time = out$survJumpC[[iEndpoint]][[iStrata]][,"time"],
                                             survival = out$survJumpC[[iEndpoint]][[iStrata]][,"survival"],
                                             dSurvival = out$survJumpC[[iEndpoint]][[iStrata]][,"dSurvival"],
                                             index_survival = out$survJumpC[[iEndpoint]][[iStrata]][,"index.survival"],
                                             index_dSurvival1 = out$survJumpC[[iEndpoint]][[iStrata]][,"index.dsurvival1"],
                                             index_dSurvival2 = out$survJumpC[[iEndpoint]][[iStrata]][,"index.dsurvival2"],
                                             lastSurv = out$lastSurv[[iEndpoint]][iStrata,2],
                                             lastdSurv = out$lastSurv[[iEndpoint]][iStrata,1], 
                                             iidNuisance = iidNuisance,
                                             nJump = NROW(out$survJumpC[[iEndpoint]][[iStrata]]))

            ls.intT <- calcIntegralSurv2_cpp(time = out$survJumpT[[iEndpoint]][[iStrata]][,"time"],
                                             survival = out$survJumpT[[iEndpoint]][[iStrata]][,"survival"],
                                             dSurvival = out$survJumpT[[iEndpoint]][[iStrata]][,"dSurvival"],
                                             index_survival = out$survJumpT[[iEndpoint]][[iStrata]][,"index.survival"],
                                             index_dSurvival1 = out$survJumpT[[iEndpoint]][[iStrata]][,"index.dsurvival1"],
                                             index_dSurvival2 = out$survJumpT[[iEndpoint]][[iStrata]][,"index.dsurvival2"],
                                             lastSurv = out$lastSurv[[iEndpoint]][iStrata,1],
                                             lastdSurv = out$lastSurv[[iEndpoint]][iStrata,2], 
                                             iidNuisance = iidNuisance,
                                             nJump = NROW(out$survJumpT[[iEndpoint]][[iStrata]]))

            ## evaluate compute integral just before the observation time (possibly shifted by tau)
            index.dSurvivalT.tau <- prodlim::sindex(jump.times = ls.intT$time,
                                                    eval.times = out$survTimeC[[iEndpoint]][[iStrata]][,"time"] - threshold[iEndpoint]) + 1
            index.dSurvivalC.0 <- prodlim::sindex(jump.times = ls.intC$time,
                                                  eval.times = out$survTimeC[[iEndpoint]][[iStrata]][,"time"]) + 1
            index.dSurvivalC.tau <- prodlim::sindex(jump.times = ls.intC$time,
                                                    eval.times = out$survTimeT[[iEndpoint]][[iStrata]][,"time"] - threshold[iEndpoint]) + 1
            index.dSurvivalT.0 <- prodlim::sindex(jump.times = ls.intT$time,
                                                  eval.times = out$survTimeT[[iEndpoint]][[iStrata]][,"time"]) + 1
            
            ## get survivals
            out$survTimeC[[iEndpoint]][[iStrata]] <- cbind(out$survTimeC[[iEndpoint]][[iStrata]],
                                                           "int.dSurvivalT-threshold_lower" = ls.intT$intSurv_lower[index.dSurvivalT.tau+1],
                                                           "int.dSurvivalT-threshold_upper" = ls.intT$intSurv_upper[index.dSurvivalT.tau+1],
                                                           "int.dSurvivalC_0_lower" = ls.intC$intSurv_lower[index.dSurvivalC.0+1],
                                                           "int.dSurvivalC_0_upper" = ls.intC$intSurv_upper[index.dSurvivalC.0+1])

            out$survTimeT[[iEndpoint]][[iStrata]] <- cbind(out$survTimeT[[iEndpoint]][[iStrata]],
                                                           "int.dSurvivalC-threshold_lower" = ls.intC$intSurv_lower[index.dSurvivalC.tau+1],
                                                           "int.dSurvivalC-threshold_upper" = ls.intC$intSurv_upper[index.dSurvivalC.tau+1],
                                                           "int.dSurvivalT_0_lower" = ls.intT$intSurv_lower[index.dSurvivalT.0+1],
                                                           "int.dSurvivalT_0_upper" = ls.intT$intSurv_upper[index.dSurvivalT.0+1])


            if(iidNuisance){
                colnames(ls.intC$intSurv_deriv) <- c("index.jump","time","index.param.surv","value.surv","index.param.dsurv1","value.dsurv1","index.param.dsurv2","value.dsurv2")
                colnames(ls.intT$intSurv_deriv) <- c("index.jump","time","index.param.surv","value.surv","index.param.dsurv1","value.dsurv1","index.param.dsurv2","value.dsurv2")
                out$survJumpC[[iEndpoint]][[iStrata]] <- ls.intC$intSurv_deriv
                out$survJumpT[[iEndpoint]][[iStrata]] <- ls.intT$intSurv_deriv
                
                out$survTimeC[[iEndpoint]][[iStrata]] <- cbind(out$survTimeC[[iEndpoint]][[iStrata]],
                                                               "index.int.dSurvivalT-threshold" = index.dSurvivalT.tau,
                                                               "indexMax.int.dSurvivalT-threshold" = NROW(ls.intT$intSurv_deriv)-1,
                                                               "index.int.dSurvivalC_0" = index.dSurvivalC.0,
                                                               "indexMax.int.dSurvivalC_0" = NROW(ls.intC$intSurv_deriv)-1)
                
                out$survTimeT[[iEndpoint]][[iStrata]] <- cbind(out$survTimeT[[iEndpoint]][[iStrata]],
                                                               "index.int.dSurvivalC-threshold" = index.dSurvivalC.tau,
                                                               "indexMax.int.dSurvivalC-threshold" = NROW(ls.intC$intSurv_deriv)-1,
                                                               "index.int.dSurvivalT_0" = index.dSurvivalT.0,
                                                               "indexMax.int.dSurvivalT_0" = NROW(ls.intT$intSurv_deriv)-1)
                
            }else{
                out$survJumpC[[iEndpoint]][[iStrata]] <- matrix(0, nrow = 0, ncol = 0)
                out$survJumpT[[iEndpoint]][[iStrata]] <- matrix(0, nrow = 0, ncol = 0)
            }
        }
    }

    ## ** export
    return(out)
    
}


## * .sindex2
## e.g. jump.times = 1:3, eval.times = c(0,1,1.1,2,3,4) should give c(1,2,2,3,4,4)
##'  .sindex2(jump.times = 1:3, eval.times = c(0,1,1.1,2,3,4))
##'  prodlim::sindex(jump.times = 1:3, eval.times = c(0,1,1.1,2,3,4))+1
##'  3 - prodlim::sindex(jump.times = 1:3, eval.times = c(0,1,1.1,2,3,4), strict = TRUE, comp = "greater") + 1     
.sindex2 <- function(jump.times, eval.times){
    return(length(jump.times)-prodlim::sindex(jump.times = jump.times, eval.times = eval.times, strict = TRUE, comp = "greater")+1)
}
        
######################################################################
### BuyseTest-Peron.R ends here
