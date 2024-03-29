% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/valid.R
\name{validFCTs}
\alias{validFCTs}
\alias{validCharacter}
\alias{validClass}
\alias{validDimension}
\alias{validInteger}
\alias{validLogical}
\alias{validNames}
\alias{validNumeric}
\alias{validPath}
\title{Check Arguments of a function.}
\usage{
validCharacter(
  value1,
  name1 = as.character(substitute(value1)),
  valid.length,
  valid.values = "character",
  refuse.NULL = TRUE,
  refuse.duplicates = FALSE,
  method = NULL,
  addPP = TRUE
)

validClass(
  value1,
  name1 = as.character(substitute(value1)),
  valid.class,
  type = "inherits",
  method = NULL,
  addPP = TRUE
)

validDimension(
  value1,
  value2 = NULL,
  name1 = as.character(substitute(value1)),
  name2 = as.character(substitute(value2)),
  valid.dimension = NULL,
  type = c("NROW", "NCOL"),
  method = NULL,
  addPP = TRUE
)

validInteger(
  value1,
  name1 = as.character(substitute(value1)),
  valid.length,
  min = NULL,
  max = NULL,
  refuse.NA = TRUE,
  refuse.NULL = TRUE,
  refuse.duplicates = FALSE,
  method = NULL,
  addPP = TRUE
)

validLogical(
  value1,
  name1 = as.character(substitute(value1)),
  valid.length,
  refuse.NULL = TRUE,
  refuse.NA = TRUE,
  method = NULL,
  addPP = TRUE
)

validNames(
  value1,
  name1 = as.character(substitute(value1)),
  refuse.NULL = TRUE,
  valid.length = NULL,
  valid.values = NULL,
  required.values = NULL,
  refuse.values = NULL,
  method = NULL,
  addPP = TRUE
)

validNumeric(
  value1,
  name1 = as.character(substitute(value1)),
  valid.length,
  valid.values = NULL,
  min = NULL,
  max = NULL,
  refuse.NA = TRUE,
  refuse.NULL = TRUE,
  refuse.duplicates = FALSE,
  method = NULL,
  addPP = TRUE,
  unlist = FALSE
)

validPath(
  value1,
  name1 = as.character(substitute(value1)),
  type,
  method = NULL,
  addPP = TRUE,
  extension = NULL,
  check.fsep = FALSE
)
}
\arguments{
\item{value1}{the value of the (first) argument to be checked}

\item{name1}{the name of the (first) argument.}

\item{valid.length}{the acceptable length(s) for the argument. If \code{NULL} no test is performed.}

\item{valid.values}{the acceptable value(s) for the argument. If \code{NULL} no test is performed. Can also be "character" or "character_or_logical".}

\item{refuse.NULL}{should an error be output if value is \code{NULL}.}

\item{refuse.duplicates}{should an error be output if value contains duplicated values.}

\item{method}{the name of the function using the argument.}

\item{addPP}{add ": " after the name of the function in the error message.}

\item{type}{For \code{validDimension}: the type of operator used to check the dimensions. For \code{validPath} either "dir" or "file" to check whether to path points to an existing directory or file.}

\item{value2}{the second value of a second argument whose dimensions should be consistent with the first one}

\item{name2}{the name of the second argument.}

\item{min}{the minimum acceptable value}

\item{max}{the maximum acceptable value}

\item{refuse.NA}{should an error be output if value contains \code{NA}.}

\item{required.values}{values that must appear in the argument}

\item{refuse.values}{values that must not appear in the argument}

\item{unlist}{[logical] flatten argument before check.}

\item{extension}{filter the files by the type of extension.}

\item{check.fsep}{display a warning when the separator is not correctly specified in}

\item{validClass}{the acceptable classes(s) for the argument.}

\item{validDimension}{the acceptable dimension for the argument. If \code{NULL} then name2 is used as a reference.}
}
\value{
An invisible \code{TRUE} or an error message.
}
\description{
Check the validity of the arguments in functions.
}
\concept{check}
\keyword{internal}
