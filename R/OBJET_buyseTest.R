#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% objet BuyseRes
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClass(
  
  Class="BuyseRes",
  
  representation(
    delta = "matrix", 
    count_favorable = "matrix",      
    count_unfavorable = "matrix",
    count_neutral = "matrix",    
    count_uninf = "matrix",
    index_neutralT="vector",
    index_neutralC = "vector",
    index_uninfT="vector",
    index_uninfC="vector",
    n_pairs = "numeric",
    delta_boot = "array", 
    p.value = "vector",    
    Delta_quantile = "matrix",
    endpoints = "vector"
  ),
  
  validity = function(object){
    #cat("--- BuyseRes : checking --- ")
    
    n.strata <- nrow(object@delta)
    n.outcome <- ncol(object@delta)
    
    if(nrow(object@count_favorable)!=n.strata || ncol(object@count_favorable)!=n.outcome)
    {stop("validity[BuyseRes] : \'@count_favorable\' does not match  \'@delta\' dimensions \n",
          "dim(@delta) : ",n.strata," ",n.outcome," \n",
          "dim(@count_favorable) : ",row(object@count_favorable)," ",ncol(object@count_favorable)," \n")
    }
    
    if(nrow(object@count_unfavorable)!=n.strata || ncol(object@count_unfavorable)!=n.outcome)
    {stop("validity[BuyseRes] : \'@count_unfavorable\' does not match  \'@delta\' dimensions \n",
          "dim(@delta) : ",n.strata," ",n.outcome," \n",
          "dim(@count_favorable) : ",row(object@count_unfavorable)," ",ncol(object@count_unfavorable)," \n")
    }
    
    if(nrow(object@count_neutral)!=n.strata || ncol(object@count_neutral)!=n.outcome)
    {stop("validity[BuyseRes] : \'@count_neutral\' does not match  \'@delta\' dimensions \n",
          "dim(@delta) : ",n.strata," ",n.outcome," \n",
          "dim(@count_neutral) : ",row(object@count_neutral)," ",ncol(object@count_neutral)," \n")
    }
    
    if(nrow(object@count_uninf)!=n.strata || ncol(object@count_uninf)!=n.outcome)
    {stop("validity[BuyseRes] : \'@count_uninf\' does not match  \'@delta\' dimensions \n",
          "dim(@delta) : ",n.strata," ",n.outcome," \n",
          "dim(@count_uninf) : ",row(object@count_uninf)," ",ncol(object@count_uninf)," \n")
    }
    
    if(length(object@index_neutralT)!=length(object@index_neutralC))
    {stop("validity[BuyseRes] : \'@index_neutralT\' does not match  \'@index_neutralC\' dimensions \n",
          "length(@index_neutralT) : ",length(object@index_neutralT)," \n",
          "length(@index_neutralC) : ",length(object@index_neutralC)," \n")
    }
    
    if(length(object@index_uninfT)!=length(object@index_uninfC))
    {stop("validity[BuyseRes] : \'@index_uninfT\' does not match  \'@index_uninfT\' dimensions \n",
          "length(@index_uninfT) : ",length(object@index_uninfT)," \n",
          "length(@index_uninfC) : ",length(object@index_uninfC)," \n")
    }
    
    #     if(object@n_pairs %% 1 != 0)
    #     {stop("validity[BuyseRes] : wrong specification of \'@n_pairs\' \n",
    #           "\'n_pairs\' must be an integer \n",
    #           "@n_pairs : ",object@n_pairs," \n")
    #     }
    
    if(length(object@p.value)!=n.outcome)
    {stop("validity[BuyseRes] : wrong specification of \'@n.outcome\' \n",
          "must have length n.outcome : ",n.outcome," \n",
          "length(@p.value) : ",length(object@p.value)," \n")
    }
    
    if(dim(object@delta_boot)[1]!=n.strata || dim(object@delta_boot)[2]!=n.outcome)
    {stop("validity[BuyseRes] : wrong specification of \'@delta_boot\' \n",
          "must have dim[1]=",n.strata," and dim[2]=",n.outcome," \n",
          "dim(object@delta_boot) : ",paste(dim(object@delta_boot),collapse=" ")," \n")
    }
    
    if(nrow(object@Delta_quantile)!=2 || ncol(object@Delta_quantile)!=n.outcome)
    {stop("validity[BuyseRes] : wrong specification of \'@Delta_quantile\' \n",
          "must have dimensions : 2 ",n.outcome," \n",
          "dim(@Delta_quantile) : ",nrow(object@Delta_quantile)," ",ncol(object@Delta_quantile)," \n")
    }
    
    if(length(object@endpoints)!=n.outcome)
    {stop("validity[BuyseRes] : wrong specification of \'@endpoints\' \n",
          "must have length n.outcome : ",n.outcome," \n",
          "length(@endpoints) : ",length(object@endpoints)," \n")
    }
    
    #cat(" : valid BuyseRes  \n")
    return(TRUE)} 
  
  
)



setMethod(f ="summary",
          signature ="BuyseRes",
          definition = function(object,type="pc",digit=3,trace=TRUE){
            
            # preparation
            if(type %in% c("nb","pc") == FALSE){
              stop("summary[BuyseRes] : wrong specification of \'type\' \n",
                   "valid types : \"nb\" \"pc\" \n",
                   "proposed \'type\' : ",type,"\n")
            }
            
            # mise en forme
            n.endpoint <- ncol(object@delta)
            n.strata <- nrow(object@delta)
            
            table <- list(nb=data.frame(matrix(NA,nrow=(n.strata+1)*n.endpoint,ncol=14)),
                          pc=data.frame(matrix(NA,nrow=(n.strata+1)*n.endpoint,ncol=14)))
            names(table$nb) <- c("endpoint_threshold","strata","n.total","n.favorable","n.unfavorable","n.neutral","n.uninf","delta","Delta","CIinf.Delta","CIsup.Delta","n.bootstrap","p.value","")
            names(table$pc) <- c("endpoint_threshold","strata","pc.total","pc.favorable","pc.unfavorable","pc.neutral","pc.uninf","delta","Delta","CIinf.Delta","CIsup.Delta","n.bootstrap","p.value","")
            Delta_boot <- apply(apply(object@delta_boot,c(2,3),sum),2,cumsum)
            if(n.endpoint==1){Delta_boot <- rbind(Delta_boot)}
            
            # global
            index.global <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1
            table$nb[index.global,"endpoint_threshold"] <- object@endpoints
            table$nb[index.global,"strata"] <- "global"
            table$nb[index.global,"n.favorable"] <- colSums(object@count_favorable)
            table$nb[index.global,"n.unfavorable"] <- colSums(object@count_unfavorable)
            table$nb[index.global,"n.neutral"] <- colSums(object@count_neutral)
            table$nb[index.global,"n.uninf"] <- colSums(object@count_uninf)
            table$nb[index.global,"delta"] <- colSums(object@delta)
            table$nb[index.global,"Delta"] <- cumsum(colSums(object@delta))
            table$nb[index.global,"CIinf.Delta"] <- cumsum(colSums(object@delta))+object@Delta_quantile["2.5%",]
            table$nb[index.global,"CIsup.Delta"] <- cumsum(colSums(object@delta))+object@Delta_quantile["97.5%",]
            table$nb[index.global,"n.bootstrap"] <- apply(Delta_boot,1,function(x){sum(!is.na(x))})
            table$nb[index.global,"p.value"] <- object@p.value
            table$nb[index.global,ncol(table$nb)] <- sapply(object@p.value,function(x){
              if(is.na(x)){NA}else if(x<0.001){"***"}else if(x<0.01){"**"}else if(x<0.05){"*"}else if(x<0.1){"."}else{""}
            })
            table$nb[index.global,"n.total"] <- rowSums(table$nb[index.global,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])
            
            table$pc[index.global,"endpoint_threshold"] <- object@endpoints
            table$pc[index.global,"strata"] <- "global"
            table$pc[index.global,"pc.favorable"] <- 100*colSums(object@count_favorable)/table$nb[1,"n.total"]#table$nb[index.global,"n.total"]
            table$pc[index.global,"pc.unfavorable"] <- 100*colSums(object@count_unfavorable)/table$nb[1,"n.total"]#table$nb[index.global,"n.total"]
            table$pc[index.global,"pc.neutral"] <- 100*colSums(object@count_neutral)/table$nb[1,"n.total"]#table$nb[index.global,"n.total"]
            table$pc[index.global,"pc.uninf"] <- 100*colSums(object@count_uninf)/table$nb[1,"n.total"]#table$nb[index.global,"n.total"]
            table$pc[index.global,"delta"] <- colSums(object@delta)
            table$pc[index.global,"Delta"] <- cumsum(colSums(object@delta))
            table$pc[index.global,"CIinf.Delta"] <- cumsum(colSums(object@delta))+object@Delta_quantile["2.5%",]
            table$pc[index.global,"CIsup.Delta"] <- cumsum(colSums(object@delta))+object@Delta_quantile["97.5%",]  
            table$pc[index.global,"p.value"] <- object@p.value
            table$pc[index.global,"n.bootstrap"] <- apply(Delta_boot,1,function(x){sum(!is.na(x))})
            table$pc[index.global,ncol(table$pc)] <- sapply(object@p.value,function(x){
              if(is.na(x)){NA}else if(x<0.001){"***"}else if(x<0.01){"**"}else if(x<0.05){"*"}else if(x<0.1){"."}else{""}
            })
            table$pc[index.global,"pc.total"] <- rowSums(table$pc[index.global,c("pc.favorable","pc.unfavorable","pc.neutral","pc.uninf")]) 
            
            for(iter_strata in 1:n.strata){
              index.strata <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1+iter_strata
              table$nb[index.strata,"endpoint_threshold"] <- object@endpoints
              table$nb[index.strata,"strata"] <- iter_strata
              table$nb[index.strata,"n.favorable"] <- object@count_favorable[iter_strata,]
              table$nb[index.strata,"n.unfavorable"] <- object@count_unfavorable[iter_strata,]
              table$nb[index.strata,"n.neutral"] <- object@count_neutral[iter_strata,]
              table$nb[index.strata,"n.uninf"] <- object@count_uninf[iter_strata,]
              table$nb[index.strata,"delta"] <- object@delta[iter_strata,]
              table$nb[index.strata,"n.total"] <- rowSums(table$nb[index.strata,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])
              
              table$pc[index.strata,"endpoint_threshold"] <- object@endpoints
              table$pc[index.strata,"strata"] <- iter_strata
              table$pc[index.strata,"pc.favorable"] <- 100*object@count_favorable[iter_strata,]/table$nb[1,"n.total"]#table$nb[index.strata,"n.total"]
              table$pc[index.strata,"pc.unfavorable"] <- 100*object@count_unfavorable[iter_strata,]/table$nb[1,"n.total"]#table$nb[index.strata,"n.total"]
              table$pc[index.strata,"pc.neutral"] <- 100*object@count_neutral[iter_strata,]/table$nb[1,"n.total"]#table$nb[index.strata,"n.total"]
              table$pc[index.strata,"pc.uninf"] <- 100*object@count_uninf[iter_strata,]/table$nb[1,"n.total"]#table$nb[index.strata,"n.total"]
              table$pc[index.strata,"delta"] <- object@delta[iter_strata,]      
              table$pc[index.strata,"pc.total"] <- rowSums(table$pc[index.strata,c("pc.favorable","pc.unfavorable","pc.neutral","pc.uninf")])
            }
            
            
            # affichage
            if(trace==TRUE){
              
              if(all(na.omit(table$nb$n.bootstrap==0))){
                indexCol <- 1:9
              }else{
                indexCol <- 1:14
              }
             
              # rounding      
              if(type=="nb"){
                table.print <- table$nb[table$nb$strata=="global",,drop=FALSE]
                param.signif <- c("delta","Delta","CIinf.Delta","CIsup.Delta")
                table.print[,param.signif] <- sapply(table.print[,param.signif],signif,digit=digit)
              }
              if(type=="pc"){
                table.print <- table$pc[table$pc$strata=="global",,drop=FALSE]
                param.round <- c("pc.total","pc.favorable","pc.unfavorable","pc.neutral","pc.uninf")
                table.print[,param.round] <- sapply(table.print[,param.round],round,digit=2)
                param.signif <- c("delta","Delta","CIinf.Delta","CIsup.Delta")
                table.print[,param.signif] <- sapply(table.print[,param.signif],signif,digit=digit)                
              }
              rownames(table.print) <- 1:nrow(table.print)
              print(table.print[,indexCol,drop=FALSE])         
            }
            
            # export
            return(invisible(table))
          }
)
