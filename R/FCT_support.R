#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% User functions
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constStrata <- function(data,strata,sep=".",lex.order = FALSE,trace=TRUE,as.numeric=FALSE){
  
  if(any(strata %in% names(data) == FALSE)){
    stop("constStrata : wrong specification of \'strata\' \n",
         "some columns requested are missing in data \n",
         "missing strata : ",paste(strata[strata %in% names(data) == FALSE],collapse=" "),"\n",
         "available variables in data : ",paste(names(data)[names(data) %in% strata == FALSE],collapse=" "),"\n")
  }
  
  resInteractions <- interaction(as.list(data[,strata]),drop = TRUE,lex.order=lex.order,sep=sep)
  levels <- levels(resInteractions)
  n.levels <- length(levels)
  
  
  # display
  if(trace==TRUE){
    table_tempo <- as.numeric(table(resInteractions))
    max.num <- 5 #nchar(max(n.levels))
    ncharLevels <- nchar(levels)
    
    textLevels <- sapply(1:n.levels,function(x){
      paste(levels[x],paste(rep(" ",max(6-ncharLevels[x],max(ncharLevels)-ncharLevels[x])),collapse="")," : ",table_tempo[x],sep="")
    })
    
    cat(n.levels," strata were founded on the ",length(strata)," strata variable",if(length(strata)>1){"s"}," (",paste(strata,collapse=" "),")\n",        
        "(",rep("#",max.num),") strata ",paste(rep(" ",max(0,max(ncharLevels)-6)),collapse=""),": number of observations \n",sep="")
    
    for(iter_level in 1:n.levels){            
      cat("(",iter_level,")",paste(rep(" ",max.num-nchar(iter_level),collapse=""))," ",textLevels[[iter_level]],"\n",sep="")            
    }
    
    cat("(total) ",rep(" ",max(ncharLevels,6))," : ",length(resInteractions),"\n",sep="")
  }
  
  # conversion
  if(as.numeric==TRUE){
    resInteractions <- as.numeric(resInteractions)
  }
  
  # export 
  return(resInteractions)
}
