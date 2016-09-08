'.onAttach' <- function(libname, pkgname="BuyseTest") {
  desc <- utils::packageDescription(pkgname)
  packageStartupMessage(desc$Package, " version ",desc$Version)
}