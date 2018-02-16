require(zoo)
require(TMB)
require(toal)


data.test <- read.HTI.RAT('~/Dropbox/PCAD4Cod/HTI data for James/AT16S012840914 (G8 4,0,48).RAT', fs = 12000)
data.test.sub <- subset(data.test, subset = seconds < 60)
head(data.test)

TagFreq.detect(detections = data.test.sub$seconds, n.tags = 4 ,
               frequencies = seq(0.9,1.1,0.001), plot = T, minGap = 0.01)

data("Pen_SeaBass")

## Speed of sound (with simualted error)
c <- 1500
c.real <- rnorm(n = 1, mean = 1500,sd = 1)

## Get track time-or-arrival values for first fish
fish.ids <- unique(Pen_SeaBass[[1]]$Period)
toa.expected <- xyz2toa(hydrophone.positions = Pen_SeaBass[[2]][,2:4],
                        tag.positions = subset(Pen_SeaBass[[1]],
                                               Period == fish.ids[1])[1:100,c(3:5,1)],
                        c = c.real)

## Simulate soruces of error
hydrophoneMovement.sd <- 0.02
pingDetection <- 0.8
hist(rnorm(n = 10000,mean = 0, sd = hydrophoneMovement.sd)*100,
     main = 'Hydrophone movement', col = 'grey', xlab = 'Centimeters')
toa.real <- toaSimError(toa = toa.expected, reciever.fs = 12000,
                        hydrophoneMovement.sd = hydrophoneMovement.sd,
                        p.detect = pingDetection, c = c.real)

# Plot differences between expected and real toa
hist((toa.expected - toa.real)*1000, main = 'Latency Error/Detection',
     xlab = 'Milliseconds', col = 'grey')

hist(apply(X = toa.real, MARGIN = 1,
           FUN = function(x){
             length(which(!is.na(x)))
           }),breaks = c(0:4),
     col = 'grey', main = 'Hydrophone Detections/Ping',xlab = '# Hydrophones')

### Testing ###
hydrophone.pos <- Pen_SeaBass$hydrophones[,2:4]
toa <- toa.real

## Remove NAs
toa.clean <- toa
toa.clean[is.na(toa)] <- -9999

data <- list(
  H = as.matrix(hydrophone.pos),
  toa = as.matrix(toa.clean),
  nh = nrow(hydrophone.pos),
  np = nrow(toa)
)

## Interlop mising times
top.interp <- spline(x = 1:data$np, y = toa[,4], xout = 1:data$np)[[2]]

params <- list(
  XYZ = as.matrix(cbind(x = runif(n = data$np, # Random position at time of ping
                                  min = min(hydrophone.pos$x),
                                  max = max(hydrophone.pos$x)),
                        y = runif(n = data$np,
                                  min = min(hydrophone.pos$y),
                                  max = max(hydrophone.pos$y)),
                        z = runif(n = data$np,
                                  min = min(hydrophone.pos$z),
                                  max = max(hydrophone.pos$z)))),
  top = top.interp,	# Estimated random times of pings
  dl = matrix(data = rnorm(n = data$np*data$nh, mean = 0,sd = 1),
              nrow = data$np,
              ncol = data$nh),		# latency error
  c  = c, # Speed of sound

  logSigma_dl = rep(-3,length = data$nh),		# Sigma for latency error

  logD_xy = -2,    		# Log SD of XY movement/unit time
  logD_z = -2,        # Log SD of Z movement/unit time

  logSigma_toa = -8, # Time of arrival SD
  logScale_toa = -3,		# scale-parameter for t-dist

  log_t_part = -3		# t-part of mixture model
)

# gdbsource('testing.r', T)

# Make optimization function
obj <- MakeADFun(data = data,
                 parameters = params,
                 random = c('XYZ', 'top', 'dl'),
                 DLL="yaps3D_Pelagic", # Position, speed of sound, time of pings
                 inner.control = list(maxit = 1000),
                 silent=F)
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

#Extract data
sd_xy <- matrix(plsd$XYZ, ncol=3)
yapsRes <- data.frame(x=pl$XYZ[,1], y=pl$XYZ[,2], z=pl$XYZ[,3],
                      top=pl$top, sd_x=sd_xy[,1], sd_y=sd_xy[,2], sd_y=sd_xy[,3])

