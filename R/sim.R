### sim.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Sep  1 2024 (18:53) 
## Version: 
## Last-Updated: jul  3 2026 (14:45) 
##           By: Brice Ozenne
##     Update #: 6
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
### sim.R ends here
