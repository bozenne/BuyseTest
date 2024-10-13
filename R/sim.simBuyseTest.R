### sim.simBuyseTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Sep  1 2024 (18:53) 
## Version: 
## Last-Updated: Sep  1 2024 (18:55) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @export
sim.simBuyseTest <- function(x, ...){

    out <- rbind(lava::sim(x$lvm[["T"]], n = x$n["C"], latent = x$latent),
                 lava::sim(x$lvm[["C"]], n = x$n["C"], latent = x$latent))
    if(!is.null(x$prefix.cluster)){
        out[[x$name.cluster]] <- paste0(x$prefix.cluster,out[[x$name.cluster]])
    }

    ## ** export
    return(out)

}

##----------------------------------------------------------------------
### sim.simBuyseTest.R ends here
