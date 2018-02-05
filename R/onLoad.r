.onLoad <- function(libname, pkgname) {
  ## Compile c++ code
  require(TMB)
  compile("./yaps.cpp")

  #### Compile and run TMB-model
  # Load compiled library
  dyn.load(dynlib("./yaps"))
}
