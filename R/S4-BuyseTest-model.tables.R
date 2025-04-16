## * model.tables (documentation)
#' @docType methods
#' @name S4BuyseTest-model.tables
#' @title Extract Summary for Class "S4BuyseTest"
#' @aliases model.tables,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Extract a summary of the results from the \code{\link{BuyseTest}} function.
#' 
#' @param x output of \code{\link{BuyseTest}}
#' @param percentage [logical] Should the percentage of pairs of each type be displayed ? Otherwise the number of pairs is displayed.
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016),
#' \code{"favorable"} displays the proportion in favor of the treatment (also called Mann-Whitney parameter), as described in Fay et al. (2018).
#' \code{"unfavorable"} displays the proportion in favor of the control.
#' Default value read from \code{BuyseTest.options()}.
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param strata [logical] should the strata-specific results be displayed or the results pooled across strata?
#' Can also be \code{NULL} to display both.
#' @param columns [character vector] subset of columns to be output (e.g. \code{"endpoint"}, \code{"favorable"}, ...).
#' Can also be \code{"summary"} or \code{"print"} to only select columns displayed in the summary or print. \code{NULL} will select all columns.
#' @param ... arguments to be passed to \code{\link{S4BuyseTest-confint}}
#'
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{S4BuyseTest-class}} for a presentation of the \code{S4BuyseTest} object. \cr
#'   \code{\link{S4BuyseTest-confint}} to output confidence interval and p-values in a matrix format.
#' 
#' @examples
#' library(data.table)
#' 
#' dt <- simBuyseTest(1e2, n.strata = 3)
#' 
#'  \dontrun{
#'  BT <- BuyseTest(treatment ~ TTE(eventtime, status = status) + Bin(toxicity), data=dt)
#'  }
#'  \dontshow{
#'  BT <- BuyseTest(treatment ~ TTE(eventtime, status = status) + Bin(toxicity), data=dt, n.resampling = 10, trace = 0)
#'  }
#'  model.tables(BT)
#'  model.tables(BT, percentage = FALSE)
#'  model.tables(BT, statistic = "winRatio")
#'
#' @keywords methods
#' @author Brice Ozenne

## * mode.tables (code)
#' @rdname S4BuyseTest-model.tables
#' @exportMethod model.tables
setMethod(f = "model.tables",
          signature = "S4BuyseTest",
          definition = function(x, percentage = TRUE, statistic = NULL, conf.level = NULL, strata = NULL,
                                columns = "summary", ...){

              ## ** normalize and check arguments
              option <- BuyseTest.options()
              mycall <- match.call()
              if(is.null(statistic)){
                  statistic <- option$statistic
              }
              
              validLogical(percentage,
                           name1 = "percentage",
                           valid.length = 1,
                           refuse.NA = FALSE, 
                           method = "model.tables[S4BuyseTest]")

              statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                  "netbenefit" = "netBenefit",
                                  "winratio" = "winRatio",
                                  "favorable" = "favorable",
                                  "unfavorable" = "unfavorable",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio","favorable","unfavorable"),
                             valid.length = 1,
                             method = "model.tables[S4BuyseTest]")

              scoring.rule <- slot(x,"scoring.rule")
              strata.level <- x@level.strata
              if(is.null(strata)){
                  if(length(strata.level)==1 || (attr(scoring.rule,"test.paired") && "strata" %in% names(mycall) == FALSE)){
                      strata <- "global"
                  }else{
                      strata <- c("global",strata.level)
                  }
              }else if(identical(strata,FALSE)){
                  strata <- "global"
              }else if(identical(strata,TRUE)){
                  strata <- strata.level
              }else if(is.numeric(strata)){
                  validInteger(strata,
                               name1 = "strata",
                               valid.length = NULL,
                               min = 1,
                               max = length(strata.level),
                               refuse.NULL = TRUE,
                               refuse.duplicates = TRUE,
                               method = "model.tables[S4BuyseTest]")
                  strata <- strata.level[strata]
              }else{
                  validCharacter(strata,
                                 name1 = "strata",
                                 valid.length = NULL,
                                 valid.values = c("global",strata.level),
                                 refuse.NULL = FALSE,
                                 method = "model.tables[S4BuyseTest]")
              }

              ## ** load info from object
              hierarchical <- slot(x,"hierarchical")
              endpoint <- slot(x,"endpoint")
              n.endpoint <- length(endpoint)
              level.strata <- slot(x,"level.strata")
              n.strata <- length(level.strata)
              n.pairs <- slot(x,"n.pairs")
              count.favorable <- slot(x,"count.favorable")
              count.unfavorable <- slot(x,"count.unfavorable")
              count.neutral <- slot(x,"count.neutral")
              count.uninf <- slot(x,"count.uninf")

              delta <- coef(x, statistic = statistic, cumulative = FALSE, strata = c("global",strata.level), simplify = FALSE)
              Delta <- coef(x, statistic = statistic, cumulative = TRUE, strata = "global", simplify = FALSE)
              n.resampling <- dim(x@nResampling)[1]

              method.inference <- x@method.inference

              ## ** compute confidence intervals and p-values
              outConfint  <- confint(x, statistic = statistic, cumulative = TRUE, strata = "global", conf.level = conf.level,  ...)

              ## ** generate summary table
              ## *** prepare
              table <- data.frame(matrix(NA,nrow=(n.strata+1)*n.endpoint,ncol=20),
                                  stringsAsFactors = FALSE)
              names(table) <- c("endpoint","restriction","threshold","weight","strata",
                                "total","favorable","unfavorable","neutral","uninf",
                                "delta","Delta","Delta(%)","information(%)",
                                "lower.ci","upper.ci","null","p.value","significance","n.resampling")
              endpoint.restriction.threshold <- rep(NA, (n.strata+1)*n.endpoint)

              ## *** global statistic
              index.global <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1
              
              table[index.global,"favorable"] <- as.double(colSums(count.favorable))
              table[index.global,"unfavorable"] <- as.double(colSums(count.unfavorable))
              table[index.global,"neutral"] <- as.double(colSums(count.neutral))
              table[index.global,"uninf"] <- as.double(colSums(count.uninf))
              table[index.global,"total"] <- rowSums(table[index.global,c("favorable","unfavorable","neutral","uninf")])

              table[index.global,"restriction"] <- x@restriction
              table[index.global,"endpoint"] <- x@endpoint
              endpoint.restriction.threshold[index.global] <- names(x@endpoint)
              table[index.global,"threshold"] <- x@threshold
              table[index.global,"weight"] <- x@weightEndpoint
              table[index.global,"strata"] <- "global"
              
              table[index.global,"delta"] <- delta["global",]
              table[index.global,"Delta"] <- Delta["global",]
              table[index.global,"Delta(%)"] <- 100*Delta["global",]/Delta["global",n.endpoint]

              ## *** strata-sepcific statistic
              for(iStrata in 1:n.strata){
                  index.strata <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1+iStrata
                  table[index.strata,"favorable"] <- count.favorable[iStrata,]
                  table[index.strata,"unfavorable"] <- count.unfavorable[iStrata,]
                  table[index.strata,"neutral"] <- count.neutral[iStrata,]
                  table[index.strata,"uninf"] <- count.uninf[iStrata,]
                  
                  table[index.strata,"strata"] <- level.strata[iStrata]
                  table[index.strata,"endpoint"] <- x@endpoint
                  endpoint.restriction.threshold[index.strata] <- names(x@endpoint)
                  table[index.strata,"threshold"] <- x@threshold
                  table[index.strata,"restriction"] <- x@restriction
                  table[index.strata,"weight"] <- x@weightEndpoint
                  table[index.strata,"delta"] <- delta[iStrata+1,]
              }

              ## *** total 
              table[,"total"] <- rowSums(table[,c("favorable","unfavorable","neutral","uninf")])

              ## *** percentage
              if(identical(percentage, TRUE)){
                  nTot.pairs <- sum(n.pairs)
              
                  table$total <- 100*table$total/nTot.pairs
                  table$favorable <- 100*table$favorable/nTot.pairs
                  table$unfavorable <- 100*table$unfavorable/nTot.pairs
                  table$neutral <- 100*table$neutral/nTot.pairs
                  table$uninf <- 100*table$uninf/nTot.pairs
              }

              ## *** information fraction and co
              table[index.global,"information(%)"] <- 100*cumsum(colSums(count.favorable+count.unfavorable)/sum(count.favorable+count.unfavorable))
             
              ## *** compute CI and p-value
              if("lower.ci" %in% names(outConfint)){
                  table[index.global,"lower.ci"] <- outConfint[,"lower.ci"]
              }
              if("upper.ci" %in% names(outConfint)){
                  table[index.global,"upper.ci"] <- outConfint[,"upper.ci"]
              }

              table[index.global,"null"] <- outConfint[,"null"]
              table[index.global,"p.value"] <- outConfint[,"p.value"]
              table[index.global,"n.resampling"] <- attr(outConfint,"n.resampling")

              ## *** restrict to strata
              if(!is.null(strata)){
                  endpoint.restriction.threshold <- endpoint.restriction.threshold[table$strata %in% strata]
                  table <- table[table$strata %in% strata,,drop = FALSE]                      
              }

              ## ** subset columns
              if(identical(columns,"summary") || any(sapply(1:9, function(iNum){identical(columns,paste0("summary",iNum))})) || identical(columns,"print")){
                  if(identical(columns,"summary")){
                      columns <- option$summary.display[[1]]
                  }else if(identical(columns,"print")){
                      columns <- option$print.display
                  }else{
                      columns <- option$summary.display[[which(paste0("summary",1:9)==columns)]]
                  }
                  columns <- setdiff(columns,"significance")
                  if("CI" %in% columns){
                      index.CI <- which(columns == "CI")
                      if(index.CI == 1){
                          columns <- c("lower","upper",columns[2:length(columns)])
                      }else if(index.CI == length(columns)){
                          columns <- c(columns[1:(index.CI-1)],"lower.ci","upper.ci")
                      }else{
                          columns <- c(columns[1:(index.CI-1)],"lower.ci","upper.ci",columns[(index.CI+1):length(columns)])
                      }
                  }
                  if(all(is.na(table$restriction))){
                      columns <- setdiff(columns,"restriction")
                  }
                  if(all(table$weight==1)){
                      columns <- setdiff(columns,"weight")
                  }
                  if(all(abs(table$threshold)<=1e-12)){
                      columns <- setdiff(columns, "threshold")
                  }
                  if(identical(strata, "global")){
                      columns <- setdiff(columns,"strata")
                  }
                  if("delta" %in% columns && "Delta" %in% columns && n.endpoint == 1 && (n.strata==1 || identical(strata, "global"))){
                      columns <- setdiff(columns, "delta")
                  }
              
                  if(method.inference == "none"){
                      columns <- setdiff(columns,c("lower.ci","upper.ci","p.value","n.resampling"))
                  }else if(attr(method.inference,"ustatistic")){
                      columns <- setdiff(columns,"n.resampling")
                  }else if(attr(method.inference,"permutation")){
                      columns <- setdiff(columns,c("lower.ci","upper.ci"))
                  }              
              }
              if(!is.null(columns)){
                  validCharacter(columns,
                                 name1 = "columns",
                                 valid.length = NULL,
                                 valid.values = names(table),
                                 refuse.NULL = FALSE,
                                 method = "summary[S4BuyseTest]")
                  table <- table[,columns,drop=FALSE]
              }

              ## ** export
              attr(table,"endpoint") <- endpoint.restriction.threshold
              attr(table,"transform") <- attr(outConfint,"nametransform")
              attr(table,"n.resampling") <- attr(outConfint,"n.resampling")
              attr(table,"method.ci.resampling") <- attr(outConfint,"method.ci.resampling")
              attr(table,"warning") <- attr(outConfint,"warning")
              return(table)
          }
)



