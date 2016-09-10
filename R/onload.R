'.onAttach' <- function(libname, pkgname="BuyseTest") {
  desc <- utils::packageDescription(pkgname)
  packageStartupMessage(desc$Package, " version ",desc$Version)
  BuyseTest.options(reinitialise = TRUE) # generate .BuyseTest-options when loading the package
}