#########################################################################
## All functions on this file are modified from the origional
## YAPS github repo: https://github.com/baktoft/yaps_sciRep.git on the date:
## 2018-02-01.  These functions are licensed under the GNU Affero General
## Public License V.3
##########################################################################

yaps2D <- function(hydrophone.pos, toa){
  ################################################################################################################################################################################################
  #	https://github.com/baktoft/yaps_sciRep
  #	Baktoft, Gjelland, Økland and Thygesen (2017). Positioning of aquatic animals based on Time-of-arrival and random walk models using YAPS (Yet Another Positioning Solver).
  #	Example code to reproduce YAPS-results presented in figure 5 from the manuscript. Position estimation is based on raw time-of-arrival data obtained from a tow-track conducted in the field.
  #	TMB needs to be downloaded and installed beforehand - see https://github.com/kaskr/adcomp/wiki/Download
  ################################################################################################################################################################################################
  #	This example file should be able to run trouble free, once TMB is installed.
  #	Tested on Windows 7 using R x64 3.3.1; Windows 10 using R x64 3.4.1; Ubuntu 16.04.3 LTS using R 3.4.1
  #	In case of problems, please contact: hba@aqua.dtu.dk or submit as an issue on github. Thanks!
  ################################################################################################################################################################################################

  #load required packages - install if necessary
  require(zoo)
  require(TMB)
  # To test TMB installation
  # runExample(all=TRUE)

  #Compile TMB-model - only needed once
  compile("./src/yaps.cpp")


  #### Compile and run TMB-model
  # Load compiled library
  dyn.load(dynlib("./src/yaps"))
  # Convert to R function
  obj <- MakeADFun(data = datTmb,
                   parameters = params,
                   DLL="yaps",
                   random=c("XY","top"), # Position, speed of sound, time of pings
                   inner.control = list(maxit = 500000),
                   silent=F)
  # Run model
  ## use iterative methods to find optimal parameters
  opt <- nlminb(start = inits,
                objective = obj$fn,
                gradient = obj$gr)

  #Obtain parameter estimates and standard deviations
  obj$fn()
  pl <- obj$env$parList()
  jointrep <- sdreport(obj, getJointPrecision=TRUE)
  param_names <- rownames(summary(jointrep))
  sds <- summary(jointrep)[,2]
  summ <- data.frame(param=param_names, sd=sds)
  plsd <- split(summ[,2], f=summ$param)

  #Extract data
  sd_xy <- matrix(plsd$XY, ncol=2)
  yapsRes <- data.frame(x=pl$XY[,1], y=pl$XY[,2], top=pl$top+T0, sd_x=sd_xy[,1], sd_y=sd_xy[,2])


  #Plot results and compare to gps and umap

  plotRes(yapsRes)
}
