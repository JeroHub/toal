#########################################################################
## Author: James Campbell, Leiden University
## Created on: 2018-03-01
## License: GPL3
##########################################################################

#' Group oxygen analysis for swimtunnel experiments
#'
#' @param oxygen A data.frame with columns: trial, group, speed, measure, temperature, and oxygenConsumption.
#' @param groupID A Trials x Group x Member array holding the index values for individuals.
#' @param params Optional parameter holding starting values for parameter estimates.
#'
#' @return A table of estimated positions in 3D.
#' @export
groupSwim <- function(oxygen, groupID, params = NULL){

  require(zoo)
  require(TMB)

  ## Transform data.frames into multi-dimentional arrays

  ## Remove NAs
  oxygen[is.na(OxygenGroup)] <- -999
  groupID[is.na(groupID)] <- -999

  ## Prepare input data
  data <- list(
    ID = as.matrix(hydrophone.pos),
    OxygenGroup = as.matrix(toa.clean),
    Temperature = nrow(hydrophone.pos),
    Speed = nrow(toa),
    nTrials = c,
    nGroups = ,
    nMembers = ,
    nSpeeds = ,
    nMeasures =
  )


  ## Prepare starting values for params
  if (length(params) != 6){
    message('Using default start parameters')
    params <- list(
      XYZ = as.matrix(start.xyz),
      top = top.interp,	# Estimated random times of pings
      dl = matrix(data = rnorm(n = data$np*data$nh, mean = 0,sd = -10),
                  nrow = data$np,
                  ncol = data$nh),		# latency error
      #c  = c, # Speed of sound
      logSigma_bi = -7.4, #log transformed SD of pulse intervals
      # logSigma_dl = rep(-10,length = data$nh),		# Sigma for latency error
      logSigma_dl = -10,		# Sigma for latency error
      logSigma_toa = -15, # Time of arrival SD
      logSigma_x = -3.5,
      logSigma_y = -3.5,
      logSigma_z = -3.5)
  }else{
    message('Using user specified start parameters')
    params <- append(
      list(XYZ = as.matrix(start.xyz),
           top = top.interp,	# Estimated random times of pings
           dl = matrix(data = rnorm(n = data$np*data$nh, mean = 0,
                                    sd = exp(params$logSigma_dl)),
                       nrow = data$np,
                       ncol = data$nh)#,		# latency error)
         #  c = c
           ),
      params)} # User specificed start parameters


  message('Running YAPS...')
  # Make optimization function
  obj <- MakeADFun(data = data,
                   parameters = params,
                   random = c('XYZ', 'top', 'dl'),
                   DLL="yaps3D_Pelagic", # Position, speed of sound, time of pings
                   inner.control = list(maxit = max.iterations),
                   silent=F, checkParameterOrder = T)
  # Run model
  ## use iterative methods to find optimal parameters
  opt <- nlminb(start = unlist(params[4:length(params)]),
                objective = obj$fn,
                gradient = obj$gr)

  # Extract results
  obj$fn()
  pl <- obj$env$parList()
  jointrep <- sdreport(obj, getJointPrecision=TRUE)
  param_names <- rownames(summary(jointrep))
  sds <- summary(jointrep)[,2]
  summ <- data.frame(param=param_names, sd=sds)
  plsd <- split(summ[,2], f=summ$param)


  return(pl)
}
