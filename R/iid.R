### iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  7 2019 (11:20) 
## Version: 
## Last-Updated: jan  7 2019 (15:18) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

iid.BuyseRes <- function(x, endpoint = NULL, ...){
    x.iid <- attr(x@covariance,"iid")
    
    if(is.null(endpoint)){
        ## iid decomposition over all endpoints
        return(apply(x.iid, MARGIN = c(1,3), sum))
    }else{
        ## iid decomposition for each endpoint
        if(any(endpoint %in% x@endpoint == FALSE)){
            txt <- endpoint[endpoint %in% x@endpoint == FALSE]
            stop("Endpoint(s) \"",paste0(txt,collapse = "\" \""),"\" not found \n")
        }else{
            return(attr(x@covariance,"iid")[endpoint])
        }
    }
    
}

######################################################################
### iid.R ends here
