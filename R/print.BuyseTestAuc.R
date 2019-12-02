### print.BuyseTestAuc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2019 (18:08) 
## Version: 
## Last-Updated: dec  2 2019 (19:19) 
##           By: Brice Ozenne
##     Update #: 15
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

print.BuyseTestAuc <- function(x, ...){
    if(NROW(x) < attr(x,"n.fold")){
        print.data.frame(x)
    }else{
        label.upper <- paste0(attr(x,"contrast")[2],">",attr(x,"contrast")[1])
        label.lower <- paste0(attr(x,"contrast")[1],">",attr(x,"contrast")[2])
        x$direction <- sapply(x$direction, function(iD){
            if(iD==">"){return(label.upper)}else if(iD=="<"){return(label.lower)}else{return(iD)}
        })
        print.data.frame(x[x$fold == "global",c("direction","estimate","se","lower","upper","p.value")], row.names = FALSE)
    }
}

######################################################################
### print.BuyseTestAuc.R ends here
