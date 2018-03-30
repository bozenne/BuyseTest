BuyseTest.env <- new.env() # create a specific environment for the package

'.onAttach' <- function(libname, pkgname="BuyseTest") {
    desc <- utils::packageDescription(pkgname)
    packageStartupMessage(desc$Package, " version ",desc$Version)
  
   BuyseTest.option(reinitialise = TRUE) # generate .BuyseTest-option when loading the package   
}
