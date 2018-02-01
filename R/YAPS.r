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
  rm(list=ls())
  set.seed(42) #for reproducible results

  #load required packages - install if necessary
  require(zoo)
  require(TMB)
  # To test TMB installation
  # runExample(all=TRUE)


  #Position of hydrophones. H14 malfunctioned and is not used - (x,y ) = (-9999,-9999)
  Hx = c(0,143.34,105.93,83.45,-86.3,-133.42,-107.29,15.05,36.91,-91.02,-44.82,
         61.99,113.24, -9999, 86.9,-125.92,21.67,52.63,55.33,-2.12,-41.71,-22.99,
         -78.85,-83.19,74.96,27.66,44.05,-115.85)
  Hy = c(0,-78.85,-85.87,-121.16,22.39,-55.03,31.69,-24.81,-34.47,4.96,6.59,
         -15.86,-53.76,-9999, -40.88,-33.9,-70.53,-93.06,-46.82,-24.19,-27.04,
         -60.77,-58.98,-19.39,-66.14,-5.73,-10.58,10.58)
  message('Hydrophones: ',length(Hx))

  #### Prepare TOAs
  # Must be in a matrix, where rows are hydrophones, and columns are pings
  toa <- read.table('./yaps_sciRep/toa.txt') #Time-of-arrival obtained from the hydrophones.
  toa <- t(as.matrix(toa))

  #### Subset TOAs
  #We suggest using a sub-sample for first testing (e.g. first 500 pings) as computation time is faster
  #first toa-observation to include
  nstart <- 1
  #number of pings to include - increase to 3244 for complete track. Be carefull if going below ~100 as number of detecting hydrophones is quite low in the beginning of the track (often <3 hydrophones detecting each ping); too many of these can cause numerical problems.
  n <- 1000
  nmax <- min(n + nstart - 1, ncol(toa))	#last toa-observation to include

  toa <- toa[,nstart:nmax]
  T0 <- min(toa, na.rm=TRUE)
  toa <- toa - T0


  #Increase pNA (range 0.00 - 1.00 ) to randomly increase probability of missed detections
  pNA <- 0.00
  toa[which(rbinom(length(toa),1,pNA) == 1)] <- NA

  #Subsample the hydrophone array by excluding hydros specified in noHs (range 1 - 28).
  noHs <- c(14) #Hydrophone 14 malfunctioned and is not used

  #This configuration leaves eight hydrophones to cover the entire study site
  # noHs <- c(14,28,16,10,5,8,26,27,19,18,13,3,20,9,25,24,21,22, 15,1)
  toa[noHs, ] <- NA


  toa[is.na(toa)] <- -9999 			#NA's are translated to -9999

  #Get data prepared for TMB
  datTmb <- list(H = matrix(c(Hx,Hy), ncol = 2), # x, y hydrophone positiosn
                 nh = length(Hx), # Number of hydrophones
                 np = ncol(toa), # Number of pings
                 toa = toa)


  ## Get list of parameters to be estimated by TMB-model
  params <- 	list(
    XY = matrix(c(mean(datTmb$H[,1]) + rnorm(ncol(datTmb$toa), sd=50), mean(datTmb$H[,2]) + rnorm(ncol(datTmb$toa), sd=50)), ncol=2),	#positions
    top = na.approx(apply(datTmb$toa, 2, function(k) {median(k[k != -9999])}), rule=2),								#time of ping
    ss=rnorm(ncol(datTmb$toa), 1415, 5),																			#speed of sound
    logD_xy = -2,				#diffusivity of transmitter movement (D_xy in ms)
    # logD_b = -15,				#diffusivity of burst interval (D_b in ms)
    logSigma_bi = -5,			#sigma  burst interval (sigma_bi in ms)
    logD_v = -3,				#diffusivity of speed of sound (D_v in ms)
    logSigma_toa = -8,			#sigma for Gaussian
    logScale = -3,				#scale parameter for t-distribution
    log_t_part = -3				#Mixture ratio between Gaussian and t
  )

  ## Get list of initial values for estimated fixed parameters
  inits <- c(-2, -5, -10, -5, -1, -2)
  names(inits) <- c('init_logD_xy', 'init_logSigma_bi', 'init_logD_v',
                    'init_logSigma_toa', 'init_logScale', 'init_log_t_part')
  inits


  #Compile TMB-model - only needed once
  compile("./src/yaps.cpp")


  #### Compile and run TMB-model
  # Load compiled library
  dyn.load(dynlib("./src/yaps"))
  # Convert to R function
  obj <- MakeADFun(data = datTmb,
                   parameters = params,
                   DLL="yaps",
                   random=c("XY","ss","top"), # Position, speed of sound, time of pings
                   inner.control = list(maxit = 500000),
                   silent=F)
  # Run model
  ## use iterative methods to find optimal parameters
  opt <- nlminb(inits,obj$fn,obj$gr)

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
