#' @name validFCTs
#' @aliases validClass
#' @aliases validDimension
#' @aliases validInteger
#' @aliases validLogical
#' @aliases validNames
#' @aliases validNumeric
#' @aliases validPath
#' @title Check argument in functions
#' 
#' @description Check the validity of the arguments in functions
#' 
#' @param value1 the value of the (first) argument to be checked
#' @param value2 the second value of a second argument whose dimensions should be consistent with the first one
#' @param name1 the name of the (first) argument.
#' @param name2 the name of the second argument.
#' @param validClass the acceptable classes(s) for the argument. 
#' @param validDimension the acceptable dimension for the argument. If \code{NULL} then name2 is used as a reference.
#' @param validLength the acceptable length(s) for the argument. If \code{NULL} no test is performed.
#' @param validValues the acceptable value(s) for the argument. If \code{NULL} no test is performed. Can also be "character" or "character_or_logical".
#' @param superClasses uses the \code{is} function instead of \code{class} to test the class of the object.
#' @param refuse.NULL should an error be output if value is \code{NULL}.
#' @param refuse.NA should an error be output if value contains \code{NA}.
#' @param refuse.duplicates should an error be output if value contains duplicated values.
#' @param refuse.values values that must not appear in the argument
#' @param type For \code{validDimension}: the type of operator used to check the dimensions. For \code{validPath} either "dir" or "file" to check whether to path points to an existing directory or file.
#' @param requiredValues values that must appear in the argument
#' @param min the minimum acceptable value
#' @param max the maximum acceptable value
#' @param method the name of the function using the argument.
#' @param addPP add ": " after the name of the function in the error message.
#' @param extension check that the file has the specified extension (excluding \code{.}).
#' @param checkFsep check whether the path ends with a separator.
#' 
#' @return An invisible \code{TRUE} or an error message.
#' 
#' @keywords function check

#' @rdname validFCTs
validCharacter <- function(value1, name1 = as.character(substitute(value1)), validLength, 
                           validValues = "character", refuse.NULL = TRUE, refuse.duplicates = FALSE, 
                           method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(is.null(value1)){
    
    if(refuse.NULL == TRUE){
      stop(method, "\'", name1, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check size
    n.value1 <- length(value1)
    
    if(!is.null(validLength) && n.value1 %in% validLength == FALSE){
      stop(method, "\'", name1, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name1, ") : ", n.value1, "\n")
    }
    
    #### check duplicates
    if(refuse.duplicates == TRUE && any(duplicated(value1))){
      stop(method, "\'", name1, "\' contains duplicated values: ", "\"",paste(unique(value1[duplicated(value1)]), collapse = "\" \""), "\" \n")
    }
    
    #### check values
    if(identical(validValues,"character")){
      
      if(any(is.character(value1) == FALSE)){
        stop(method, "\'", name1, "\' must be a ", if(n.value1 == 1){"character"}else{"vector of characters"}," \n", 
             "is(", name1, ") : ", paste(is(value1), collapse = " "), "\n")
      }
      
    } else if(identical(validValues,"character_or_logical")){
      
      if(any( (is.character(value1) == FALSE) * (is.logical(value1) == FALSE) > 0 )){
        stop(method, "\'", name1, "\' must be a ", if(n.value1 == 1){"character or logical"}else{"vector of characters or logicals"}," \n", 
             "is(", name1, ") : ", paste(is(value1), collapse = " "), "\n")
      }
      
    } else if(!is.null(validValues) && any(value1 %in% validValues == FALSE)){
      
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "valid values for \'", name1, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value1 %in% validValues == FALSE)>1){"s"}," for \'", name1, "\' : \"", paste(value1[value1 %in% validValues == FALSE], collapse = "\" \""), "\"\n")
      
    }
    
  }
  
  return(invisible(TRUE))
  
}

#' @rdname validFCTs
validClass <- function(value1, name1 = as.character(substitute(value1)), validClass, 
                       superClasses = TRUE, method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(superClasses == TRUE){
    
    if( all(is(value1) %in% validClass == FALSE) ){
      stop(method, "class of \'", name1, "\' must be one of the following \"", paste(validClass,collapse="\" \""), "\"  \n", 
           "proposed superclass : \"", paste(is(value1),collapse="\" \""), "\" \n")
    }  
    
  }else{
    
    if( class(value1) %in% validClass == FALSE){
      stop(method, "class of \'", name1, "\' must be \"", paste(validClass,collapse="\" \""),"\"  \n", 
           "proposed class : ", class(value1)[[1]], "\n")
    }  
    
  }
  
  return(invisible(TRUE))
  
}

#' @rdname validFCTs
validDimension <- function(value1, value2 = NULL, name1 = as.character(substitute(value1)), name2 = as.character(substitute(value2)),
                           validDimension = NULL,
                           type = c("NROW","NCOL"), method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  n.type <- length(type)
  
  #### dimension 1
  testDimension <- sapply(1:n.type, function(x){
    do.call(type[x], list(value1))
  })
  
  #### dimension 2
  
  
  if(is.null(validDimension)){
    
    validDimension <- sapply(1:n.type, function(x){
      do.call(type[x], list(value2))
    })
    test.validDimension <- TRUE
    
  }else if(is.null(name2)){
    
    test.validDimension <- FALSE
    
  }else{
    
    test.validDimension <- TRUE
    
  }
  
  #### main
  for(iter_type in 1:n.type){
    
    if(testDimension[iter_type] != validDimension[iter_type]){
      
      if(test.validDimension){
        stop(method, "dimension mismatch between argument \'", name1, "\' and argument \'", name2, "\' \n", 
             type[iter_type],"(", name1, ") = ", testDimension[iter_type], " \n", 
             type[iter_type],"(", name2, ") = ", validDimension[iter_type], " \n")  
      }else{
        stop(method, "dimension mismatch between argument \'", name1, "\' and argument \'", name2, "\' \n", 
             type[iter_type],"(", name1, ") = ", testDimension[iter_type], " \n", 
             type[iter_type],"(", name2, ") = ", validDimension[iter_type], " \n")
        
      }
      
    }
    
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
validInteger <- function(value1, name1 = as.character(substitute(value1)), validLength, 
                         validValues = NULL, min = NULL, max = NULL, 
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, 
                         method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  validNumeric(value1 = value1, name1 = name1, validLength = validLength, min = min, max = max, 
               validValues = validValues, refuse.NA = refuse.NA, refuse.NULL = refuse.NULL, refuse.duplicates = refuse.duplicates, 
               method = method)
  
  #### check integer
  if(!is.null(value1) && any(value1 %% 1 > 0)){
    stop(method, "\'", name1, "\' must contain integers not doubles \n",        
         "invalid value(s) in ", name1, " : ", paste(value1[value1 %% 1 > 0], collapse = " "), "\n")
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
validLogical <- function(value1, name1 = as.character(substitute(value1)), validLength, 
                         refuse.NULL = TRUE, refuse.NA = TRUE, 
                         method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(is.null(value1)){
    
    #### NULL
    if(refuse.NULL == TRUE){
      stop(method, "\'", name1, "\' must be logical ",if(refuse.NA == FALSE){"or NA"}," and not NULL \n")
    }
    
  }else{ 
    
    #### Size
    if(!is.null(validLength) && length(value1) %in% validLength == FALSE){
      stop(method, "\'", name1, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name1, ") : ", length(value1), "\n")
    } 
    
    #### Type
    if(any(is.logical(value1) == FALSE)){
      stop(method, "\'", name1, "\' must be ", if(refuse.NULL == FALSE){"NULL or "}, if(refuse.NA == FALSE){"NA or "},"TRUE or FALSE \n",        
           "is(", name1, ") : ", paste(is(value1), collapse = " "), "\n")
    }
    
    if(refuse.NA == TRUE && any(is.na(value1)) ){
      stop(method, "\'", name1, "\' must be logical ",if(refuse.NULL == FALSE){"or NULL"}," and not NA \n")
    }
    
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
validNames <- function(value1, name1 = as.character(substitute(value1)), refuse.NULL = TRUE,
                       validLength = NULL, validValues = NULL, requiredValues = NULL, refuse.values = NULL,
                       method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  ## type
  if(is.matrix(value1)){
    value1 <- colnames(value1)
  }
  
  if(data.table::is.data.table(value1) || is.data.frame(value1) || is.list(value1)){
    value1 <- names(value1)
  }
  
  ## tests
  if(is.null(value1)){
    
    if(refuse.NULL == TRUE){
      stop(method, "names of \'", name1, "\' must not be NULL \n")
    }
    
  }else{
    #### check size
    n.value1 <- length(value1)
    
    if(!is.null(validLength) && n.value1 %in% validLength == FALSE){
      stop(method, "\'", name1, "\' must have ", paste(validLength, collapse = " or ")," names  \n", 
           "length(names(", name1, ")) : ", n.value1, "\n")
    }
    
    #### check content
    
    if(!is.null(requiredValues) && any(requiredValues %in% value1 == FALSE)){
      
      stop(method, "\'", name1, "\' must contains: \"",paste(requiredValues[requiredValues %in% value1 == FALSE], collapse = "\" \""),"\" \n", 
           "proposed names : \"", paste(value1, collapse = "\" \""), "\"\n")  
      
    }
    
    if(!is.null(validValues) && any(value1 %in% validValues == FALSE)){
      
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "valid names for \'", name1, "\' : \"",paste(validValues, collapse = "\" \""),"\" \n", 
           "refused names : \"", paste(value1[value1 %in% validValues == FALSE], collapse = " "), "\"\n")  
      
    }
    
    if(!is.null(refuse.values) && any(value1 %in% refuse.values)){
      
      stop(method, "\'", name1, "\' contains forbidden names:", paste(value1[value1 %in% refuse.values], collapse = " "), "\"\n")  
      
    }
    
    if(any(duplicated(value1))){
      stop(method, "\'", name1, "\' must not contain duplicated names \n", 
           "duplicated names : \"", paste(value1[duplicated(value1)], collapse = " "), "\"\n")  
    }
    
  }
  
  return(invisible(TRUE))
  
}
#' @rdname validFCTs
validNumeric <- function(value1, name1 = as.character(substitute(value1)), validLength,
                         validValues = NULL , min = NULL, max = NULL,
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, 
                         method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(is.null(value1)){
    
    if(refuse.NULL == TRUE){
      stop(method, "\'", name1, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check length
    if(!is.null(validLength) && length(value1) %in% validLength == FALSE){
      stop(method, "\'", name1, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name1, ") : ", length(value1), "\n")
    }
    
    #### check NA
    if(refuse.NA == TRUE && any(is.na(value1))){
      stop(method, "\'", name1, "\' must not contain any NA \n", 
           "index of NA values : ", paste(which(is.na(value1)), collapse = " "), "\n")
    }
    
    #### check numeric
    if(any( (is.numeric(value1) == FALSE) * (is.na(value1) == FALSE) )){
      stop(method, "\'", name1, "\' must be numeric \n",        
           "is(", name1, ") : ", paste(is(value1), collapse = " "), "\n")
    }
    
    #### check duplicates
    if(refuse.duplicates == TRUE && any(duplicated(value1))){
      stop(method, "\'", name1, "\' contains duplicated values: ", paste(unique(value1[duplicated(value1)]), collapse = " "), "\n")
    }
    
    #### check min value1
    if(!is.null(min) && any(stats::na.omit(value1) < min)){
      stop(method, "\'", name1, "\' must be bigger than ", min, " \n",        
           "invalid value(s): ", paste(value1[stats::na.omit(value1) < min], collapse = " "), "\n")
    }
    
    #### check max value1
    if(!is.null(max) && any(stats::na.omit(value1) > max)){
      stop(method, "\'", name1, "\' must be smaller than ", max, " \n",        
           "invalid value(s): ", paste(value1[stats::na.omit(value1) > max], collapse = " "), "\n")
    }
    
    #### check valid values
    if(!is.null(validValues) && any(value1 %in% validValues == FALSE)){
      
      stop(method, "\'", name1, "\' contain invalid values \n", 
           "valid values for \'", name1, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value1 %in% validValues == FALSE)>1){"s"}," for \'", name1, "\' : \"", paste(value1[value1 %in% validValues == FALSE], collapse = " "), "\"\n")
      
    }
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
validPath <- function(value1, name1 = as.character(substitute(value1)), type,
                      method = NULL, addPP = TRUE, extension = NULL, checkFsep = FALSE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  validCharacter(type, validLength = 1, validValues = c("file", "dir"))
  
  try_path <- switch(type,
                     file = file.exists(value1),
                     dir = dir.exists(value1)
  )
  
  if(try_path == FALSE){
    stop(method, "\'", name1, "\' does not lead to an existing ",switch(type,"file"="file","dir"="directory")," \n", 
         "proposed value: \"", value1, "\"\n", 
         "current path: ", getwd(), "\n")
  }
  
  
  if(type == "dir"){ 
    if(checkFsep == TRUE && substr(value1, start = nchar(value1), stop = nchar(value1)) != "/"){
      warning(method, "possible bad specification of \'", name1, "\' \n", 
              "it should end with a fsep (e.g. \"/\") \n")
    }
  }else if(type == "file" && !is.null(extension)){
    fileExtension <- tools::file_ext(value1) 
    if(fileExtension %in% extension == FALSE){
      stop(method, "\'", name1, "\' has not the expected extension \n", 
           "proposed extension: \"", fileExtension, "\" \n", 
           "expected extension: \"", paste(extension, collapse = "\" \""), "\"\n")
    }
  }
  
  return(invisible(TRUE))
}