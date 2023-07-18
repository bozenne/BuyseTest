### BuyseTest-Peron.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 12 2020 (11:10) 
## Version: 
## Last-Updated: jul 18 2023 (12:03) 
##           By: Brice Ozenne
##     Update #: 548
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calcPeron
#' @noRd
calcPeron <- function(data,
                      model.tte, fitter, args,
                      method.score, paired,
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
                      level.strata,
                      n.strata,
                      strata,
                      threshold,
                      restriction,
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

        txt.fitter <- sapply(fitter, switch,
                             "prodlim" = "prodlim::Hist",
                             "survreg" = "survival::Surv",
                             NA)
        if(paired){
            txt.modelUTTE <- paste0(txt.fitter,"(",endpoint.UTTE,",",status.UTTE,") ~ ",treatment)        
        }else{
            txt.modelUTTE <- paste0(txt.fitter,"(",endpoint.UTTE,",",status.UTTE,") ~ ",treatment," + ..strata..")        
        }
        
    }else{
        tofit <- FALSE
    }
    
    ## fit survival model and prepare for extracting survival
    for(iUTTE in 1:D.UTTE){ ## iUTTE <- 1
        ls.indexAssociatedEndpoint[[iUTTE]] <-  which(endpoint == endpoint.UTTE[iUTTE])
        test.CR[iUTTE] <- any(method.score[ls.indexAssociatedEndpoint[[iUTTE]]]=="CRPeron")

        if(fitter[iUTTE]=="prodlim"){

            if(tofit){
                model.tte[[iUTTE]] <- do.call(prodlim::prodlim, args = c(list(as.formula(txt.modelUTTE[iUTTE]), data = data, discrete.level = 1e5), args))
            }

        }else if(fitter[iUTTE]=="survreg"){

            if(tofit){
                model.tte[[iUTTE]] <- do.call(survival::survreg, args = c(list(as.formula(txt.modelUTTE[iUTTE]), data = data), args))
            }
            
        }
        model.tte[[iUTTE]] <- BuyseTTEM(model.tte[[iUTTE]], treatment = treatment, level.treatment = level.treatment, level.strata = level.strata, iid = iidNuisance)
    }

    ## ** estimate quantities for scoring pairs
    for(iUTTE in 1:D.UTTE){ ## iUTTE <- 1
        iN.CR <- model.tte[[iUTTE]]$peron$n.CR
        
        ## fitted survival at event timepoints
        for(iStrata in 1:n.strata){

            for(iTreat in level.treatment){
                iStoreJump <- c("survJumpC","survJumpT")[(iTreat==level.treatment[2])+1]
                iStoreTime <- c("survTimeC","survTimeT")[(iTreat==level.treatment[2])+1]
                iStoreP <- c("p.C","p.T")[(iTreat==level.treatment[2])+1]

                iTreat.num <- as.numeric(iTreat==level.treatment[2])                
                iTime <- data[list(iTreat.num,iStrata),.SD[[endpoint.UTTE[iUTTE]]],on=c(treatment,"..strata..")]
                iTime.jump <- predict(model.tte[[iUTTE]], time = "jump", strata = iStrata, treatment = iTreat)

                iRestriction <- restriction[ls.indexAssociatedEndpoint[[iUTTE]][1]]
                if(!is.na(iRestriction)){ 
                    iTime <- pmin(iTime, iRestriction) ## take minimum between outcome and restriction
                    iTime.jump <- iTime.jump[iTime.jump <= iRestriction] ## remove jumps after restriction                     
                }
                if(length(iTime.jump)==0){iTime.jump <- 0}

                iPred1.C <- predict(model.tte[[iUTTE]], time = iTime, treatment = level.treatment[1], strata = iStrata, cause = 1)
                iPred1.T <- predict(model.tte[[iUTTE]], time = iTime, treatment = level.treatment[2], strata = iStrata, cause = 1)
                if(iN.CR>1){
                    iPred2 <- predict(model.tte[[iUTTE]], time = iTime, treatment = iTreat, strata = iStrata, cause = 2)
                }
                iPred1.iTreat.beforeJump <- predict(model.tte[[iUTTE]], time = iTime.jump-zeroPlus, treatment = iTreat,  strata = iStrata, iid = iidNuisance)
                iPred1.iTreat.afterJump <- predict(model.tte[[iUTTE]], time = iTime.jump+zeroPlus, treatment = iTreat, strata = iStrata, iid = iidNuisance) ## technically already computed in the previous lines

                iLastEstimate <- sapply(1:iN.CR, function(iCause){
                    if(is.na(iRestriction) || model.tte[[iUTTE]]$peron$last.time[iStrata,iTreat]<=iRestriction){ 
                        return(predict(model.tte[[iUTTE]], time = "last", strata = iStrata, treatment = iTreat, cause = iCause))
                    }else{ ## no remainder term if end of the survival curve after restriction (i.e. fully known survival up to the restriction)
                        return(0)
                    }
                })

                for(iEndpoint in ls.indexAssociatedEndpoint[[iUTTE]]){
                    iThreshold <- threshold[iEndpoint]

                    ## last estimate of the survival/cif
                    out$lastSurv[[iEndpoint]][iStrata,seq(from=iTreat.num+1, by = 2, length=iN.CR)] <- iLastEstimate
                    if(test.CR[iUTTE]){

                        ## *** CIF at jump times
                        iPred1.iOther.beforeTau <- predict(model.tte[[iUTTE]], time = iTime.jump - iThreshold,
                                                           treatment = setdiff(level.treatment,iTreat), strata = iStrata)
                        iPred1.iOther.afterTau <- predict(model.tte[[iUTTE]], time = iTime.jump + iThreshold,
                                                          treatment = setdiff(level.treatment,iTreat), strata = iStrata)
                        out[[iStoreJump]][[iEndpoint]][[iStrata]] <- cbind("time" = iTime.jump,
                                                                           "CIF1-threshold" = iPred1.iOther.beforeTau$cif,
                                                                           "CIF1+threshold" = iPred1.iOther.afterTau$cif,
                                                                           "dCIF" = iPred1.iTreat.afterJump$cif - iPred1.iTreat.beforeJump$cif,
                                                                           "index.CIF1-threshold" = iPred1.iOther.beforeTau$index,
                                                                           "index.CIF1+threshold" = iPred1.iOther.afterTau$index,
                                                                           "index.dCIF11" = iPred1.iTreat.beforeJump$index,
                                                                           "index.dCIF12" = iPred1.iTreat.afterJump$index)

                        if(iidNuisance){
                            out$iid[[iStoreJump]][[iUTTE]][[iStrata]] <- cbind(lava::iid(model.tte[[iUTTE]], strata = iStrata, treatment = iTreat, cause = 1),
                                                                               lava::iid(model.tte[[iUTTE]], strata = iStrata, treatment = iTreat, cause = 2))
                            out[[iStoreP]][iStrata, iEndpoint] <- NCOL(out$iid[[iStoreJump]][[iUTTE]][[iStrata]])
                        }

                        ## *** CIF at observation time (+/- threshold)
                        iPred1.C.beforeTau <- predict(model.tte[[iUTTE]], time = iTime - iThreshold,
                                                    treatment = level.treatment[1], strata = iStrata)
                        iPred1.C.afterTau <- predict(model.tte[[iUTTE]], time = iTime + iThreshold,
                                                   treatment = level.treatment[1], strata = iStrata)
                        iPred1.T.beforeTau <- predict(model.tte[[iUTTE]], time = iTime - iThreshold,
                                                    treatment = level.treatment[2], strata = iStrata)
                        iPred1.T.afterTau <- predict(model.tte[[iUTTE]], time = iTime + iThreshold,
                                                   treatment = level.treatment[2], strata = iStrata)

                        out[[iStoreTime]][[iEndpoint]][[iStrata]] <- cbind("time" = iTime, ## 0
                                                                           "CIF1C-threshold" = iPred1.C.beforeTau$cif, ## 1
                                                                           "CIF1C_0" = iPred1.C$cif, ## 2
                                                                           "CIF1C+threshold" = iPred1.C.afterTau$cif, ## 3
                                                                           "CIF1T-threshold" = iPred1.T.beforeTau$cif, ## 4
                                                                           "CIF1T_0" = iPred1.T$cif, ## 5
                                                                           "CIF1T+threshold" = iPred1.T.afterTau$cif, ## 6
                                                                           "CIF2_0" = iPred2$cif, ## 7
                                                                           "index.CIF1C-threshold" = iPred1.C.beforeTau$index, ## 8
                                                                           "index.CIF1C_0" = iPred1.C$index, ## 9
                                                                           "index.CIF1C+threshold" = iPred1.C.afterTau$index, ## 10
                                                                           "index.CIF1T-threshold" = iPred1.T.beforeTau$index, ## 11
                                                                           "index.CIF1T_0" = iPred1.T$index, ## 12
                                                                           "index.CIF1T+threshold" = iPred1.T.afterTau$index, ## 13
                                                                           "index.CIF2_0" = iPred2$index) ## 14
                    }else{

                        ## *** survival at jump time
                        iTimeTau.jump <- iTime.jump + iThreshold
                        if(!is.na(iRestriction)){ ## remove jump that are such that jump+threshold are after restriction
                            iSubset.restriction <- which(iTimeTau.jump<=iRestriction)
                        }else{
                            iSubset.restriction <- 1:length(iTimeTau.jump)
                        }

                        
                        if(iidNuisance){
                            out$iid[[iStoreJump]][[iUTTE]][[iStrata]] <- lava::iid(model.tte[[iUTTE]], strata = iStrata, treatment = iTreat)
                            out[[iStoreP]][iStrata, iEndpoint] <- NCOL(out$iid[[iStoreJump]][[iUTTE]][[iStrata]])
                            if(any(is.na(out$iid[[iStoreJump]][[iUTTE]][[iStrata]]))){ stop("NA in the iid decomposition of the survival model. \n") }
                        }

                        if(length(iSubset.restriction)==0){
                            out[[iStoreJump]][[iEndpoint]][[iStrata]] <- cbind(time = 0, ## jump time
                                                                               survival = 1, 
                                                                               dSurvival = 0,
                                                                               index.survival = NA, ## index of the survival parameter at t+\tau
                                                                               index.dsurvival1 = NA, ## index of the survival parameter before the jump
                                                                               index.dsurvival2 = NA) ## index of the survival parameter after the jump
                        }else{
                            iSurvTau.jump <- predict(model.tte[[iUTTE]], time = iTimeTau.jump[iSubset.restriction], treatment = setdiff(level.treatment, iTreat),
                                                     strata = iStrata, iid = iidNuisance)
                            out[[iStoreJump]][[iEndpoint]][[iStrata]] <- cbind(time = iTime.jump[iSubset.restriction], ## jump time
                                                                               survival = iSurvTau.jump$survival, 
                                                                               dSurvival = iPred1.iTreat.afterJump$survival[iSubset.restriction] - iPred1.iTreat.beforeJump$survival[iSubset.restriction],
                                                                               index.survival = iSurvTau.jump$index, ## index of the survival parameter at t+\tau
                                                                               index.dsurvival1 = iPred1.iTreat.beforeJump$index[iSubset.restriction], ## index of the survival parameter before the jump
                                                                               index.dsurvival2 = iPred1.iTreat.afterJump$index[iSubset.restriction]) ## index of the survival parameter after the jump
                        }

                        ## *** survival at observation time (+/- threshold)
                        iPred.C.beforeTau <- predict(model.tte[[iUTTE]], time = iTime - iThreshold, treatment = level.treatment[1], strata = iStrata)
                        iPred.C.afterTau <- predict(model.tte[[iUTTE]], time = iTime + iThreshold, treatment = level.treatment[1], strata = iStrata)

                        iPred.T.beforeTau <- predict(model.tte[[iUTTE]], time = iTime - iThreshold, treatment = level.treatment[2], strata = iStrata)
                        iPred.T.afterTau <- predict(model.tte[[iUTTE]], time = iTime + iThreshold, treatment = level.treatment[2], strata = iStrata)

                        out[[iStoreTime]][[iEndpoint]][[iStrata]] <- cbind("time" = iTime,
                                                                           "survivalC-threshold" = iPred.C.beforeTau$survival,
                                                                           "survivalC_0" = iPred1.C$survival,
                                                                           "survivalC+threshold" = iPred.C.afterTau$survival,
                                                                           "survivalT-threshold" = iPred.T.beforeTau$survival,
                                                                           "survivalT_0" = iPred1.T$survival,
                                                                           "survivalT+threshold" = iPred.T.afterTau$survival,
                                                                           "index.survivalC-threshold" = iPred.C.beforeTau$index,
                                                                           "index.survivalC_0" = iPred1.C$index,
                                                                           "index.survivalC+threshold" = iPred.C.afterTau$index,
                                                                           "index.survivalT-threshold" = iPred.T.beforeTau$index,
                                                                           "index.survivalT_0" = iPred1.T$index,
                                                                           "index.survivalT+threshold" = iPred.T.afterTau$index
                                                                           )
                    } ## End if CR
                } ## End if iEndpoint
            } ## End if treatment
        } ## End if strata
    } ## End if UTTE

    ## ** pre-compute integrals 
    for(iEndpoint in 1:length(endpoint)){ ## iEndpoint <- 1

        if(!precompute || method.score[iEndpoint] %in% c("continuous","gaussian")){next} ## only relevant for survival/ competing risk with Peron
    
        for(iStrata in 1:n.strata){  ## iStrata <- 1

            if(method.score[iEndpoint]=="SurvPeron"){
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
                ## e.g. jump.times = 1:3, eval.times = c(0,1,1.1,2,3,4) should give c(1,2,2,3,4,4)
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
                                                               "int.dSurvivalT-threshold_lower" = ls.intT$intSurv_lower[index.dSurvivalT.tau],
                                                               "int.dSurvivalT-threshold_upper" = ls.intT$intSurv_upper[index.dSurvivalT.tau],
                                                               "int.dSurvivalC_0_lower" = ls.intC$intSurv_lower[index.dSurvivalC.0],
                                                               "int.dSurvivalC_0_upper" = ls.intC$intSurv_upper[index.dSurvivalC.0])

                out$survTimeT[[iEndpoint]][[iStrata]] <- cbind(out$survTimeT[[iEndpoint]][[iStrata]],
                                                               "int.dSurvivalC-threshold_lower" = ls.intC$intSurv_lower[index.dSurvivalC.tau],
                                                               "int.dSurvivalC-threshold_upper" = ls.intC$intSurv_upper[index.dSurvivalC.tau],
                                                               "int.dSurvivalT_0_lower" = ls.intT$intSurv_lower[index.dSurvivalT.0],
                                                               "int.dSurvivalT_0_upper" = ls.intT$intSurv_upper[index.dSurvivalT.0])
            }else{
                stop("precompute terms for competing risks not implemented")
                ## to be done
            }

            if(iidNuisance){
                if(method.score[iEndpoint]=="SurvPeron"){
                    colnames(ls.intC$intSurv_deriv) <- c("index.jump","time","index.param.surv","value.surv","index.param.dsurv1","value.dsurv1","index.param.dsurv2","value.dsurv2")
                    colnames(ls.intT$intSurv_deriv) <- c("index.jump","time","index.param.surv","value.surv","index.param.dsurv1","value.dsurv1","index.param.dsurv2","value.dsurv2")
                    out$survJumpC[[iEndpoint]][[iStrata]] <- ls.intC$intSurv_deriv
                    out$survJumpT[[iEndpoint]][[iStrata]] <- ls.intT$intSurv_deriv

                    out$survTimeC[[iEndpoint]][[iStrata]] <- cbind(out$survTimeC[[iEndpoint]][[iStrata]],
                                                                   "index.int.dSurvivalT-threshold" = index.dSurvivalT.tau-1,
                                                                   "indexMax.int.dSurvivalT-threshold" = NROW(ls.intT$intSurv_deriv)-1,
                                                                   "index.int.dSurvivalC_0" = index.dSurvivalC.0-1,
                                                                   "indexMax.int.dSurvivalC_0" = NROW(ls.intC$intSurv_deriv)-1)
                
                    out$survTimeT[[iEndpoint]][[iStrata]] <- cbind(out$survTimeT[[iEndpoint]][[iStrata]],
                                                                   "index.int.dSurvivalC-threshold" = index.dSurvivalC.tau-1,
                                                                   "indexMax.int.dSurvivalC-threshold" = NROW(ls.intC$intSurv_deriv)-1,
                                                                   "index.int.dSurvivalT_0" = index.dSurvivalT.0-1,
                                                                   "indexMax.int.dSurvivalT_0" = NROW(ls.intT$intSurv_deriv)-1)
                }else{
                    ## to be done
                }
            }else{
                out$survJumpC[[iEndpoint]][[iStrata]] <- matrix(0, nrow = 0, ncol = 0)
                out$survJumpT[[iEndpoint]][[iStrata]] <- matrix(0, nrow = 0, ncol = 0)
            }
        }
    }

    ## ** export
    return(out)
    
}


## ## * .sindex2
## ## e.g. jump.times = 1:3, eval.times = c(0,1,1.1,2,3,4) should give c(1,2,2,3,4,4)
## ##'  .sindex2(jump.times = 1:3, eval.times = c(0,1,1.1,2,3,4))
## ##'  prodlim::sindex(jump.times = 1:3, eval.times = c(0,1,1.1,2,3,4))+1
## ##'  3 - prodlim::sindex(jump.times = 1:3, eval.times = c(0,1,1.1,2,3,4), strict = TRUE, comp = "greater") + 1     
## .sindex2 <- function(jump.times, eval.times){
##     return(length(jump.times)-prodlim::sindex(jump.times = jump.times, eval.times = eval.times, strict = TRUE, comp = "greater")+1)
## }
        
######################################################################
### BuyseTest-Peron.R ends here
