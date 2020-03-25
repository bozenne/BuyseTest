## * Documentation - constStrata
#' @title Strata creation
#' @name constStrata
#' 
#' @description Create strata from several variables.
#' 
#' @param data [data.frame] dataset.
#' @param strata [character vector] A vector of the variables capturing the stratification factors.
#' @param sep [character] string to construct the new level labels by joining the constituent ones.
#' @param lex.order [logical] Should the order of factor concatenation be lexically ordered ?
#' @param trace [logical] Should the execution of the function be traced ?
#' @param as.numeric [logical] Should the strata be converted from factors to numeric?
#' 
#' @details 
#' This function uses the \code{interaction} function from the \emph{base} package to form the strata.
#' 
#' @return A \emph{factor vector} or a \emph{numeric vector}.
#' 
#' @examples
#' library(data.table)
#' 
#' data(veteran,package="survival")
#'   
#' # strata with two variables : celltype and karno
#' veteran$strata1 <- constStrata(veteran,c("celltype","karno"))
#' table(veteran$strata1)
#'   
#' # strata with three variables : celltype, karno and age dichotomized at 60 years
#' veteran$age60 <- veteran$age>60
#' veteran$age60 <- factor(veteran$age60,labels=c("<=60",">60")) # convert to factor with labels
#' veteran$strata2 <- constStrata(veteran,c("celltype","karno","age60"))
#' table(veteran$strata2) # factor strata variable 
#'   
#' veteran$strata2 <- constStrata(veteran,c("celltype","karno","age60"), as.numeric=TRUE)
#' table(veteran$strata2) # numeric strata variable
#' 
#' @keywords function
#' @author Brice Ozenne

## * Function constStrata
#' @rdname constStrata
#' @export
constStrata <- function(data,strata,sep=".",lex.order = FALSE,trace=TRUE,as.numeric=FALSE){
  
  if(any(strata %in% names(data) == FALSE)){
    stop("constStrata : wrong specification of \'strata\' \n",
         "some columns requested are missing in data \n",
         "missing strata : ",paste(strata[strata %in% names(data) == FALSE],collapse=" "),"\n",
         "available variables in data : ",paste(names(data)[names(data) %in% strata == FALSE],collapse=" "),"\n")
  }
  
  if(data.table::is.data.table(data)){
    resInteractions <- data[,interaction(.SD[[1]],drop = TRUE,lex.order=lex.order,sep=sep), .SDcols = strata]
  }else{
    resInteractions <- interaction(as.list(data[,strata]),drop = TRUE,lex.order=lex.order,sep=sep)
  }
  levels <- levels(resInteractions)
  n.levels <- length(levels)
  
  
  ## ** display
  if(trace==TRUE){
    table_tempo <- as.numeric(table(resInteractions))
    max.num <- 5 #nchar(max(n.levels))
    ncharLevels <- nchar(levels)
    
    textLevels <- sapply(1:n.levels,function(x){
      paste(levels[x],paste(rep(" ",max(6-ncharLevels[x],max(ncharLevels)-ncharLevels[x])),collapse="")," : ",table_tempo[x],sep="")
    })
    
    cat(n.levels," strata were founded on the ",length(strata)," strata variable",if(length(strata)>1){"s"}," (",paste(strata,collapse=" "),")\n",        
        "(",rep("#",max.num),") strata ",paste(rep(" ",max(0,max(ncharLevels)-6)),collapse=""),": number of observations \n",sep="")
    
    for(iLevel in 1:n.levels){            
      cat("(",iLevel,")",paste(rep(" ",max.num-nchar(iLevel),collapse=""))," ",textLevels[[iLevel]],"\n",sep="")            
    }
    
    cat("(total) ",rep(" ",max(ncharLevels,6))," : ",length(resInteractions),"\n",sep="")
  }
  
  ## ** conversion
  if(as.numeric==TRUE){
    resInteractions <- as.numeric(resInteractions)
  }
  
  ## ** export 
  return(resInteractions)
}
