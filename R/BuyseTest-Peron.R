### BuyseTest-Peron.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 12 2020 (11:10) 
## Version: 
## Last-Updated: May 31 2026 (00:52) 
##           By: Brice Ozenne
##     Update #: 827
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
                      scoring.rule, args, model.tte, 
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
                      level.strata,
                      grid.strata,
                      strata,
                      threshold,
                      multiplicative.threshold,
                      restriction,
                      precompute,
                      iidNuisance,
                      out){

    zeroPlus <- 1e-12

    ## ** prepare
    ## fit model for the cumulative incidence function (survival case or competing risk case)
    if(scoring.rule %in% c("peron","efron","latta")){ 
        fitter <- "prodlim"
    }else{ ## Peron
        fitter <- "survreg"
        args$dist <- scoring.rule
    }
    efron <- as.numeric(scoring.rule == "efron")
    latta <- as.numeric(scoring.rule == "latta")

    model.tte <- fitTTEM(model.tte, data = data, BuyseTTEM = TRUE,
                         strata = strata, level.strata = level.strata, treatment = treatment, level.treatment = level.treatment,
                         endpoint.UTTE = endpoint.UTTE, status.UTTE = status.UTTE, D.UTTE = D.UTTE,
                         fitter = fitter, efron = efron, latta = latta, iidNuisance = iidNuisance, args = args)

    ## endpoint with competing risk
    ls.indexAssociatedEndpoint <- setNames(vector(mode = "list", length = D.UTTE), endpoint.UTTE)
    test.CR <- setNames(vector(mode = "logical", length = D.UTTE), endpoint.UTTE)
    for(iUTTE in 1:D.UTTE){ ## iUTTE <- 1
        ls.indexAssociatedEndpoint[[iUTTE]] <-  which(endpoint == endpoint.UTTE[iUTTE])
        test.CR[iUTTE] <- any(method.score[ls.indexAssociatedEndpoint[[iUTTE]]]=="CRPeron")
    }

    ## ** estimate quantities for scoring pairs
    n.strata <- NROW(grid.strata)
    grid.strataTreatment <- data[, .SD[!duplicated(.SD$..strata..)], by = treatment]  ## unique strata & treatment combinations


    for(iUTTE in 1:D.UTTE){ ## iUTTE <- 1
        iN.CR <- model.tte[[iUTTE]]$peron$n.CR
                
        for(iStrata in 1:n.strata){

            ## *** treatment group specific strata variable (may differ when using standardization: M,F -> M.M,M.F,F.M,F.F)
            iStrata.num <- stats::setNames(c(grid.strata[iStrata,1]+1,grid.strata[iStrata,2]+1), level.treatment)            
            ## strata variable in the time to event model
            ## may not match the strata variable in the GPC formula when the user specifies its own survival model

            if(length(model.tte[[iUTTE]]$peron$level.strata)==1){
                iStrata.model <- stats::setNames(c(1,1), level.treatment)                
            }else{
                iStrata.model <- list(as.data.frame(grid.strataTreatment[list(0,iStrata.num[1]), on = c(treatment,"..strata.."), .SD, .SDcols = c(strata,"..strata..")]),
                                      as.data.frame(grid.strataTreatment[list(1,iStrata.num[2]), on = c(treatment,"..strata.."), .SD, .SDcols = c(strata,"..strata..")]))
                names(iStrata.model) <- level.treatment
                if(any(sapply(iStrata.model,NROW)>2)){
                    stop("Cannot handle a survival model conditional on covariates that are not included in the GPC procedure as group or strata variables. \n")
                }
            }


            for(iTreat in level.treatment){

                iStoreJump <- c("survJumpC","survJumpT")[(iTreat==level.treatment[2])+1]
                iStoreTime <- c("survTimeC","survTimeT")[(iTreat==level.treatment[2])+1]
                iStoreP <- c("p.C","p.T")[(iTreat==level.treatment[2])+1]

                iTreat.num <- as.numeric(iTreat==level.treatment[2])

                ## *** Restriction
                ## contrained to be specific to a time to event variable and not to a endpoint
                iRestriction <- restriction[ls.indexAssociatedEndpoint[[iUTTE]][1]]

                ## *** Observation time
                iTime <- data[list(iTreat.num,iStrata.num[iTreat]),.SD[[endpoint.UTTE[iUTTE]]],on=c(treatment,"..strata..")]
                ## evaluate observation times beyond restriction at restriction
                if(!is.na(iRestriction)){
                    iTime <- pmin(iTime, iRestriction)
                }

                ## *** Survival/CIF at jump times
                iPred1.iTreatJump <- predict(model.tte[[iUTTE]], type = "change", treatment = iTreat,  strata = iStrata.model[[iTreat]])
                ## range(iPred1.iTreatJump$time - predict(model.tte[[iUTTE]], type = "jump", treatment = iTreat, strata = iStrata.model[[iTreat]]))
                iTime.jump <- iPred1.iTreatJump$time
                
                ## *** Survival/CIF at observed event times in both groups
                ## at observation times
                iPred1.C <- predict(model.tte[[iUTTE]], time = iTime, treatment = level.treatment[1], strata = iStrata.model[[1]], cause = 1)
                iPred1.T <- predict(model.tte[[iUTTE]], time = iTime, treatment = level.treatment[2], strata = iStrata.model[[2]], cause = 1)
                if(iN.CR>1){
                    iPred2 <- predict(model.tte[[iUTTE]], time = iTime, treatment = iTreat, strata = iStrata.model[[iTreat]], cause = 2)
                }
                                
                ## *** Survival/CIF at last observation time
                iLastEstimate <- vector(mode = "numeric", length = iN.CR)
                iLastEstimate[1] <- predict(model.tte[[iUTTE]], type = "last", strata = iStrata.model[[iTreat]], treatment = iTreat, cause = 1)
                if(iN.CR>1){
                    iLastEstimate[2] <- predict(model.tte[[iUTTE]], type = "last", strata = iStrata.model[[iTreat]], treatment = iTreat, cause = 2)
                }
                if(!is.na(iRestriction)){ ## no remainder term if end of the survival/CIF curve after restriction (i.e. fully known survival/CIF up to the restriction)
                    iLast.time <- predict(model.tte[[iUTTE]], type = "last.time", strata = iStrata.model[[iTreat]], treatment = iTreat, cause = 1)
                    iLastEstimate[1] <- (iLast.time<=iRestriction)*iLastEstimate[1]
                    if(iN.CR>1){
                        iLast.time2 <- predict(model.tte[[iUTTE]], type = "last.time", strata = iStrata.model[[iTreat]], treatment = iTreat, cause = 2)
                        iLastEstimate[2] <- (iLast.time2<=iRestriction)*iLastEstimate[2]
                    }
                }

                for(iEndpoint in ls.indexAssociatedEndpoint[[iUTTE]]){

                    ## *** apply additive or multiplicative threshold
                    iThreshold <- threshold[iEndpoint]
                    if(multiplicative.threshold[iEndpoint]){
                        iTime.jump_PLUS_threshold <- iTime.jump * iThreshold
                        iTime.jump_MINUS_threshold <- iTime.jump / iThreshold
                        iTime_PLUS_threshold <- iTime * iThreshold
                        iTime_MINUS_threshold <- iTime / iThreshold
                    }else{
                        iTime.jump_PLUS_threshold <- iTime.jump + iThreshold
                        iTime.jump_MINUS_threshold <- iTime.jump - iThreshold
                        iTime_PLUS_threshold <- iTime + iThreshold
                        iTime_MINUS_threshold <- iTime - iThreshold
                    }

                    ## *** only consider jumps (+tau) before restriction times
                    if(!is.na(iRestriction)){ 
                        iSubset.restriction <- which(iTime.jump_PLUS_threshold <= iRestriction)
                    }else{
                        iSubset.restriction <- 1:length(iTime.jump_PLUS_threshold)
                    }

                    ## *** store last estimate of the survival/cif
                    out$lastSurv[[iEndpoint]][iStrata,seq(from=iTreat.num+1, by = 2, length=iN.CR)] <- iLastEstimate
                    
                    ## competing risk vs. survival case
                    if(test.CR[iUTTE]){

                        ## *** store CIF at jump times
                        if(length(iSubset.restriction)==0){
                            out[[iStoreJump]][[iEndpoint]][[iStrata]] <- cbind("time" = 0,
                                                                               "CIF1-threshold" = 0, 
                                                                               "CIF1+threshold" = 0,
                                                                               "dCIF" = 0,
                                                                               "index.CIF1-threshold" = NA,
                                                                               "index.CIF1+threshold" = NA,
                                                                               "index.dCIF11" = NA,
                                                                               "index.dCIF12" = NA)
                        }else{
                            iPred1.iOther.beforeTau <- predict(model.tte[[iUTTE]], time = iTime.jump_MINUS_threshold[iSubset.restriction],
                                                               treatment = setdiff(level.treatment,iTreat), strata = iStrata.model[[setdiff(level.treatment,iTreat)]])
                            iPred1.iOther.afterTau <- predict(model.tte[[iUTTE]], time = iTime.jump_PLUS_threshold[iSubset.restriction],
                                                              treatment = setdiff(level.treatment,iTreat), strata = iStrata.model[[setdiff(level.treatment,iTreat)]])
                            out[[iStoreJump]][[iEndpoint]][[iStrata]] <- cbind("time" = iTime.jump[iSubset.restriction],
                                                                               "CIF1-threshold" = iPred1.iOther.beforeTau$cif,
                                                                               "CIF1+threshold" = iPred1.iOther.afterTau$cif,
                                                                               "dCIF" = iPred1.iTreatJump$cif[iSubset.restriction],
                                                                               "index.CIF1-threshold" = iPred1.iOther.beforeTau$index,
                                                                               "index.CIF1+threshold" = iPred1.iOther.afterTau$index,
                                                                               "index.dCIF11" = iPred1.iTreatJump$indexBefore[iSubset.restriction],
                                                                               "index.dCIF12" = iPred1.iTreatJump$indexBefore[iSubset.restriction])
                        }
                        
                        if(iidNuisance){
                            out$iid[[iStoreJump]][[iUTTE]][[iStrata]] <- cbind(lava::iid(model.tte[[iUTTE]], treatment = iTreat, strata = iStrata.model[[iTreat]], cause = 1),
                                                                               lava::iid(model.tte[[iUTTE]], treatment = iTreat, strata = iStrata.model[[iTreat]], cause = 2))
                            out[[iStoreP]][iStrata, iEndpoint] <- NCOL(out$iid[[iStoreJump]][[iUTTE]][[iStrata]])
                        }

                        ## *** store CIF at observation time (+/- threshold)
                        iPred1.C.beforeTau <- predict(model.tte[[iUTTE]], time = iTime_MINUS_threshold,
                                                      treatment = level.treatment[1], strata = iStrata.model[[level.treatment[1]]])
                        iPred1.C.afterTau <- predict(model.tte[[iUTTE]], time = iTime_PLUS_threshold,
                                                     treatment = level.treatment[1], strata = iStrata.model[[level.treatment[1]]])
                        iPred1.T.beforeTau <- predict(model.tte[[iUTTE]], time = iTime_MINUS_threshold,
                                                      treatment = level.treatment[2], strata = iStrata.model[[level.treatment[2]]])
                        iPred1.T.afterTau <- predict(model.tte[[iUTTE]], time = iTime_PLUS_threshold,
                                                     treatment = level.treatment[2], strata = iStrata.model[[level.treatment[2]]])

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

                        ## *** store survival at jump time
                        if(length(iSubset.restriction)==0){
                            out[[iStoreJump]][[iEndpoint]][[iStrata]] <- cbind(time = 0, ## jump time
                                                                               survival = 1, 
                                                                               dSurvival = 0,
                                                                               index.survival = NA, ## index of the survival parameter at t+\tau
                                                                               index.dsurvival1 = NA, ## index of the survival parameter before the jump
                                                                               index.dsurvival2 = NA) ## index of the survival parameter after the jump
                        }else{
                            iSurvTau.jump <- predict(model.tte[[iUTTE]], time = iTime.jump_PLUS_threshold[iSubset.restriction], treatment = setdiff(level.treatment, iTreat),
                                                     strata = iStrata.model[[setdiff(level.treatment, iTreat)]])
                            out[[iStoreJump]][[iEndpoint]][[iStrata]] <- cbind(time = iTime.jump[iSubset.restriction], ## jump time
                                                                               survival = iSurvTau.jump$survival, 
                                                                               dSurvival = iPred1.iTreatJump$survival[iSubset.restriction],
                                                                               index.survival = iSurvTau.jump$index, ## index of the survival parameter at t+\tau
                                                                               index.dsurvival1 = iPred1.iTreatJump$indexBefore[iSubset.restriction], ## index of the survival parameter before the jump
                                                                               index.dsurvival2 = iPred1.iTreatJump$indexAfter[iSubset.restriction]) ## index of the survival parameter after the jump
                        }

                        ## *** store influence function of the survival estimator
                        if(iidNuisance){
                            out$iid[[iStoreJump]][[iUTTE]][[iStrata]] <- lava::iid(model.tte[[iUTTE]], strata = iStrata.model[[iTreat]], treatment = iTreat)
                            out[[iStoreP]][iStrata, iEndpoint] <- NCOL(out$iid[[iStoreJump]][[iUTTE]][[iStrata]])
                            if(any(is.na(out$iid[[iStoreJump]][[iUTTE]][[iStrata]]))){ stop("NA in the iid decomposition of the survival model. \n") }
                        }
                        
                        ## *** survival at observation time (+/- threshold)
                        iPred.C.beforeTau <- predict(model.tte[[iUTTE]], time = iTime_MINUS_threshold, treatment = level.treatment[1], strata = iStrata.model[[level.treatment[1]]])
                        iPred.C.afterTau <- predict(model.tte[[iUTTE]], time = iTime_PLUS_threshold, treatment = level.treatment[1], strata = iStrata.model[[level.treatment[1]]])

                        iPred.T.beforeTau <- predict(model.tte[[iUTTE]], time = iTime_MINUS_threshold, treatment = level.treatment[2], strata = iStrata.model[[level.treatment[2]]])
                        iPred.T.afterTau <- predict(model.tte[[iUTTE]], time = iTime_PLUS_threshold, treatment = level.treatment[2], strata = iStrata.model[[level.treatment[2]]])

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
                if(multiplicative.threshold[iEndpoint]){
                    survTimeC_MINUS_threshold <- out$survTimeC[[iEndpoint]][[iStrata]][,"time"] / threshold[iEndpoint]
                    survTimeT_MINUS_threshold <- out$survTimeT[[iEndpoint]][[iStrata]][,"time"] / threshold[iEndpoint]
                }else{
                    survTimeC_MINUS_threshold <- out$survTimeC[[iEndpoint]][[iStrata]][,"time"] - threshold[iEndpoint]
                    survTimeT_MINUS_threshold <- out$survTimeT[[iEndpoint]][[iStrata]][,"time"] - threshold[iEndpoint]
                }
                
                index.dSurvivalT.tau <- prodlim::sindex(jump.times = ls.intT$time,
                                                        eval.times = survTimeC_MINUS_threshold) + 1
                index.dSurvivalC.0 <- prodlim::sindex(jump.times = ls.intC$time,
                                                      eval.times = out$survTimeC[[iEndpoint]][[iStrata]][,"time"]) + 1
                index.dSurvivalC.tau <- prodlim::sindex(jump.times = ls.intC$time,
                                                        eval.times = survTimeT_MINUS_threshold) + 1
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


## * modelUTTE
fitTTEM <- function(model.tte, data, BuyseTTEM,
                    strata, level.strata, treatment, level.treatment, endpoint.UTTE, status.UTTE, D.UTTE,
                    fitter, efron, latta, iidNuisance, args){

    if("n.grid" %in% names(args)){
        n.grid <- args$n.grid
        args$n.grid <- NULL
    }else{
        n.grid <- NULL
    }
    if(length(fitter)<D.UTTE & length(fitter)==1){
        fitter <- rep(fitter, D.UTTE)
    }
    
    ## ** prepare formula
    if(is.null(model.tte)){
        model.tte <- vector(length = D.UTTE, mode = "list")
        names(model.tte) <- endpoint.UTTE
        tofit <- TRUE
        if(length(args)==0){args <- NULL}

        txt.fitter <- sapply(fitter, switch,
                             "prodlim" = "prodlim::Hist",
                             "survreg" = "survival::Surv",
                             NA)

        if(is.null(strata) || attr(strata,"match")){
            if(latta){
                txt.modelUTTE <- paste0(txt.fitter,"(",endpoint.UTTE,",",status.UTTE,") ~ 1")        
            }else{
                txt.modelUTTE <- paste0(txt.fitter,"(",endpoint.UTTE,",",status.UTTE,") ~ ",treatment)        
            }            
        }else{
            if(latta){
                txt.modelUTTE <- paste0(txt.fitter,"(",endpoint.UTTE,",",status.UTTE,") ~  ..strata..")        
            }else{
                txt.modelUTTE <- paste0(txt.fitter,"(",endpoint.UTTE,",",status.UTTE,") ~ ",treatment," + ..strata..")        
            }            
        }
        if(efron & any(fitter!="prodlim")){
            stop("Can only use the Efron scoring rule with the Kaplan-Meier estimator. \n",
                 "Consider using the function \'efronlim\' to specify the survival model(s). \n")
        }
        
    }else{
        tofit <- FALSE
    }

    ## ** fit survival model and prepare for extracting survival
    for(iUTTE in 1:D.UTTE){ ## iUTTE <- 1

        if(fitter[iUTTE]=="prodlim"){

            if(tofit){
                if(efron){
                    model.tte[[iUTTE]] <- do.call(efronlim, args = c(list(as.formula(txt.modelUTTE[iUTTE]), data = data, discrete.level = 1e5), args))
                }else{
                    model.tte[[iUTTE]] <- do.call(prodlim::prodlim, args = c(list(as.formula(txt.modelUTTE[iUTTE]), data = data, discrete.level = 1e5), args))
                }
                
            }else{
                if(efron && is.null(model.tte[[iUTTE]]$efron)){
                    stop("Use function \'efronlim\' instead of \'prodlim\' when providing the survival model with argument scoring.rule=\"Efron\". \n")
                }else if(!efron && !is.null(model.tte[[iUTTE]]$efron)){
                    stop("Use function \'prodlim\' instead of \'efronlim\' when providing the survival model with argument scoring.rule=\"Peron\" or scoring.rule=\"Latta\". \n")
                }
            }

        }else if(fitter[iUTTE]=="survreg"){

            if(tofit){
                model.tte[[iUTTE]] <- do.call(survival::survreg, args = c(list(as.formula(txt.modelUTTE[iUTTE]), data = data), args))
             }
            
        }

        if(BuyseTTEM){
            model.tte[[iUTTE]] <- BuyseTTEM(model.tte[[iUTTE]], treatment = treatment,
                                            level.treatment = level.treatment,
                                            level.strata = list(NULL,level.strata)[[tofit+1]], ## only pass the original strata level when the model is fit internally
                                            iid = iidNuisance, n.grid = n.grid)
            model.tte[[iUTTE]]$efron <- efron & (fitter[iUTTE]=="prodlim")
        }
        
    }


    ## ** export
    return(model.tte)

}

######################################################################
### BuyseTest-Peron.R ends here
