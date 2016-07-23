#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% objet BuyseSample
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClass(
  
  Class="BuyseSample",
  
  representation(
    alpha="numeric",
    alternative="character",
    KMimputation="logical",
    n.bootstrap="numeric",
    type="vector",
    threshold="vector",
    power="array",
    delta = "array", 
    censoring.rate = "array", 
    p.value = "array",    
    DeltaCI ="array"
  ),
   
  validity = function(object){
    #cat("--- BuyseSample : checking --- ")
    
    D <- length(object@type)
    type <- initType(object@type,method="validity[BuyseSample]")
    D.TTE <- sum(object@type==3)
    
    #### alternative
    if(object@alternative %in% c("two.sided","less","greater") == FALSE){
      stop("validity[BuyseRes] :\'@alternative\' is not valid \n",
           "valid \'@alternative\' :  \"two.sided\" \"less\" \"greater\" \n",
           "proposed \'@alternative\' = ",object@alternative," \n")
    }
    
    #### threshold
    if(length(object@threshold)!=D){
      stop("validity[BuyseRes] :\'@threshold\' length does not match \'@type\' length  \n",
           "length(@type) = ",D," \n",
           "length(@threshold) = ",length(object@threshold)," \n")
    }
    
    #### delta
    dim_tempo <- dim(object@delta)
    
    if(length(dim_tempo)!=3){
      stop("validity[BuyseRes] : wrong specification of \'@delta\' \n",
           "\'@delta\' must have 3 dimensions \n",
           "length(dim(@delta)) = ",length(dim_tempo)," \n")
    }
    
    if(dim_tempo[2]!=D){
      stop("validity[BuyseRes] : \'@delta\' dimensions 2 does not match \'@type\' length \n",
           "length(@type) : ",type," \n",
           "dim(@delta)[2] = ",dim_tempo[2]," \n")
    }
    
    n.simul <- dim_tempo[1]
    n.n <- dim_tempo[3]
    
    names <- dimnames(object@delta)
    endpoint <- names[[2]]
    n <- names[[3]]
    
    #### power
    dim_tempo <- dim(object@power)
    
    if(length(dim_tempo)!=2){
      stop("validity[BuyseRes] : wrong specification of \'@power\' \n",
           "\'@power\' must have 2 dimensions \n",
           "length(dim(@power)) = ",length(dim_tempo)," \n")
    }
    
    if(dim_tempo[1]!=D || dim_tempo[2]!=n.n){
     stop("validity[BuyseRes] : \'@power\' does not match  \'@delta\' dimensions (2,3) \n",
          "valid dimensions : ",D," ",n.n," \n",
          "dim(@power) : ",paste(dim_tempo,collapse=" ")," \n")
    }
    
    #### censoring.rate
    dim_tempo <- dim(object@censoring.rate)
      
    if(D.TTE>0 && length(dim_tempo)!=4){
      stop("validity[BuyseRes] : wrong specification of \'@censoring.rate\' \n",
           "\'@censoring.rate\' must have 4 dimensions \n",
           "length(dim(@censoring.rate)) = ",length(dim_tempo)," \n")
    }
    
    if(D.TTE>0 &&  (dim_tempo[1]!=n.simul || dim_tempo[2]!=length(unique(endpoint[type==3])) || dim_tempo[3]!=n.n || dim_tempo[4]!=2)){
     stop("validity[BuyseRes] : wrong specification of \'@censoring.rate\' \n",
          "valid dimensions : ",n.simul," ",length(unique(endpoint[type==3]))," ",n.n," 2 \n",
          "dim(@censoring.rate) : ",paste(dim_tempo,collapse=" ")," \n")
    }
    
    #### p.value
    dim_tempo <- dim(object@p.value)
      
    if(length(dim_tempo)!=3){
      stop("validity[BuyseRes] : wrong specification of \'@p.value\' \n",
           "\'@p.value\' must have 3 dimensions \n",
           "length(dim(@p.value)) = ",length(dim(p.value))," \n")
    }
    
    if(dim_tempo[1]!=n.simul || dim_tempo[2]!=D ||dim_tempo[3]!=n.n){
      stop("validity[BuyseRes] : \'@p.value\' does not match  \'@delta\' dimensions \n",
           "valid dimensions : ",n.simul," ",D," ",n.n," \n",
           "dim(@p.value) : ",paste(dim_tempo,collapse=" ")," \n")
    }
    
    #### DeltaCI
    dim_tempo <- dim(object@DeltaCI)
    
    if(length(dim_tempo)!=3){
      stop("validity[BuyseRes] : wrong specification of \'@DeltaCI\' \n",
           "\'@DeltaCI\' must have 3 dimensions \n",
           "length(dim(@DeltaCI)) = ",length(dim(p.value))," \n")
    }
    
    if(dim_tempo[1]!=n.simul || dim_tempo[2]!=4 || dim_tempo[3]!=n.n){
      stop("validity[BuyseRes] : wrong specification of \'@DeltaCI\' \n",
           "valid dimensions : ",n.simul," 4 ",n.n," \n",
           "dim(@DeltaCI) : ",paste(dim_tempo,collapse=" ")," \n")
    }
    
    #cat(" : valid BuyseSample  \n")
    return(TRUE)} 
  
  
)



setMethod(f ="summary",
          signature ="BuyseSample",
          definition = function(object,Delta=FALSE,delta=FALSE,censoring.rate=FALSE,digit=3,trace=TRUE){
          
            
            #### preparation ####
            type <- object@type
            D <- length(type)            
            D.TTE <- sum(type==3)
            if(D.TTE==0){
              censoring.rate <- FALSE
            }
            
            dim_tempo <- dim(object@delta)
            n.simul <- dim_tempo[1]
            n.n <- dim_tempo[3]
            
            names <- dimnames(object@delta)
            endpoint <- names[[2]]
            n <- names[[3]]
            
            threshold <- object@threshold
            threshold[is.na(threshold)] <- ""
            
            #### table ####            
            table <- data.frame(matrix(NA,nrow=5+D+2*length(unique(endpoint[type==3])),ncol=n.n))
            names(table) <- n
            rownames_table <-  c("power","Delta (median)","Delta_inf (median)","Delta_sup (median)","n.bootstrap (median)",
                                 paste(endpoint,".delta.",threshold," (median)",sep="")
            )
            
            if(D.TTE>0){
              rownames_table <- c(rownames_table,
                                  paste(unique(endpoint[type==3]),".censoringT (median)",sep=""),
                                  paste(unique(endpoint[type==3]),".censoringC (median)",sep="")
                                  )
            }
                            
            rownames(table) <- rownames_table
          
            table["power",] <- object@power[D,]         
            table[c("Delta (median)","Delta_inf (median)","Delta_sup (median)","n.bootstrap (median)"),] <- apply(object@DeltaCI,c(2,3),median)
            table[paste(endpoint,".delta.",threshold," (median)",sep=""),] <- apply(object@delta,c(2,3),median)            
            if(D.TTE>0){
            table[paste(unique(endpoint[type==3]),".censoringT (median)",sep=""),] <- apply(object@censoring.rate[,,,"treatment",drop=FALSE],c(2,3),median)
            table[paste(unique(endpoint[type==3]),".censoringC (median)",sep=""),] <- apply(object@censoring.rate[,,,"control",drop=FALSE],c(2,3),median)
            }
            
          #### display ####

          if(trace==TRUE){
            index_table <- "power"
            if(Delta==TRUE){
              index_table <- c(index_table,"Delta (median)","Delta_inf (median)","Delta_sup (median)","n.bootstrap (median)")
            }
            if(delta==TRUE){
              index_table <- c(index_table,paste(endpoint,".delta.",threshold," (median)",sep=""))
            }
            if(censoring.rate==TRUE){
              index_table <- c(index_table,
                               paste(unique(endpoint[type==3]),".censoringT (median)",sep=""),
                               paste(unique(endpoint[type==3]),".censoringC (median)",sep=""))
            }
            
            table_display <- table[index_table,,drop=FALSE]
            table_display <- round(table_display,digit=digit)
            names(table_display) <- paste(names(table_display)," ",sep="")
            
            cat("   # alpha = ",object@alpha," (",object@alternative,") \n",sep="")    
            cat("   # n.simul = ",n.simul," | n.bootstrap = ",object@n.bootstrap," \n",sep="")    
            cat("   # endpoints \n")
            for(iter_m in 1:D){
              cat("      >  endpoint ",iter_m," = \"",endpoint[iter_m],"\" - type = \"",c("binary","continuous","timeToEvent")[type[iter_m]],"\"",sep="")
              if(type[iter_m] %in% c(2,3)){cat(" , threshold =",threshold[iter_m])}            
              cat("\n")
            }
            if(any(type==3)){cat("   # KM imputation : ",object@KMimputation,"\n")}
            cat("   # results \nn")
            print(table_display)
          }
          
          #### export ####
          return(invisible(table))
          }
)
