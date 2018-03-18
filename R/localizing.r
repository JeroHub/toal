#' Time-of-arrival difference localization
#'
#' Localize sound sources using an array of 4 or more hydrophones.
#' Four hydrophones will result in 1 degree of ambiguity and two position estimates, while more hydrophones will result in a single estimate.
#' This is an implementation by the spherical interpolation methods described by Schau & Robinson (1987) and Smith & Abel (1987) which provides a closed form solution for estimating source locations from the relative differences between time-of-arrival at each hydrophone.
#'
#' The parameter \code{toa} must be a N column matrix or data.frame where each column represents a unique hydrophone in the array and each row represents a detection.
#' The values in this parameter represent the time of arrival in seconds at each hydrophone.
#' The parameter \code{hydrophone.positions} should be a N by 3 matrix or data.frame, where each row represents a hydrophone (in the same order as appearing in the columns of \code{toa}) and the three columns are the x, y, and z positions in meters.
#'
#' When analyzing a four hydrophone array, two results will be given by two forms of the error minimizing formula.
#' These forms are indicaed by the column `eq` in the resulting data.frame.
#' In some setups, one form (either `+` or `-`) will always result in the corret location.
#' If this is the case, you can use this column to easily filter out the erronious locations.
#'
#' @param toa data.frame or matrix holding times of arrival for each hydrophone  See details for formatting requirements.
#' @param hydrophone.positions data.frame or matrix holding information about hydrophone positions.  See details for formatting.
#' @param c Speed of sound in the medium.  The default value represents a general approximation of sea water.
#'
#' @return A data.frame holding the `id` of the localization, which is the row number from the `toa` parameter highlighting the TOA differences used to make the estimate and the resulting location of the estimate.
#' @export
#'
#' @examples
#'
#' @references Smith, J. O., Abel, J. S. 1987. Closed-Form Least-Squares Source Location Estimation from Range-Difference Measurements.  IEEE TRANSACTIONS ON ACOUSTICS, SPEECH, AND SIGNAL PROCESSING. 35(12):1661-1669.
#'
#' Schau, H. C., Robinson, A. Z. 1987. Passive Source Localization Employing Intersecting Spherical Surfaces from Time-of-Arrival Differences.  IEEE TRANSACTIONS ON ACOUSTICS, SPEECH, AND SIGNAL PROCESSING. 35(8):1223-1225
TOA.localization <- function(toa, hydrohpone.positions, c = 1500){

  N <- nrow(hydrohpone.positions)
  if(N<4){
    stop('Less than 4 rows (hydrophone positions) detected')
  }
  if(ncol(toa) != N){
    stop('Number of columns in toa match number of rows in hydrophone.positions')
  }

  ## Nested function for calculating distance between points in 3d space
  dist <- function(a,b){
    ## Calc distance between two locations
    dif <- abs(a - b)
    return(sqrt(sum(dif^2)))
  }

  ## Define S: sensor locations, use first hydrophone (row 1) as origin
  # Normalize locations so H1 is at 0,0,0
  hydrohpone.positions.norm <- apply(hydrohpone.positions,MARGIN = 1,FUN = function(x){
    x - hydrohpone.positions[1,]
  })
  hydrohpone.positions.norm <- do.call(what = 'rbind',args = hydrohpone.positions.norm)

  # x,y,z position of sensors 2-4 (sensor = receivers)
  S  = as.matrix(hydrohpone.positions.norm[2:N,])

  ## Loop through each TOA and estimate a source position
  if(N==4){
    len <- nrow(toa)*2
  }else{
    len <- nrow(toa)
  }
  estimates <- data.frame(id = vector(mode = 'integer',length = len),
                          x = vector(mode = 'double',length = len),
                          y = vector(mode = 'double',length = len),
                          z = vector(mode = 'double',length = len),
                          error = vector(mode = 'double',length = len),
                          eq = vector(mode = 'character',length = len),
                          stringsAsFactors = F)
  for(i in 1:nrow(toa)){
    ## Define d: Range distance differences between all sensors
    # TOA differences are proportional to range differences, convert TOA diffs (s) to range diffs (m)
    # Create matrix of distances between sensors i,j
    d <- matrix(nrow = N-1,ncol = 1)
    for(j in 1:(N-1)){
      TOA.diff <- toa[i,j+1] - toa[i,1]
      range.diff <- TOA.diff*c # distance differences
      d[j,1] <- range.diff
    }

    ## Define R: Distances between Reference sensor (id: 1) and other sensors
    R <- matrix(nrow = N-1,ncol = 1)
    for(j in 1:(N-1)){
      R[j,1] <- sqrt(sum(S[j,]^2))
    }

    ## delta: R - d
    delta <- R^2 - d^2

    if(N > 4){
      #### Use generalized form of localization
      ## 'Weighting matrix.
      # Set all values to 1 and assume all sesnors are equally reliable
      W = diag(x = rep.int(1,times = N-1))

      ## Identity matrix
      I = diag(x = rep.int(1,times = N-1))

      ## Orthogonal compliment of P, with respect to d
      P.orth_d <- I - ((d %*% t(d)) / as.vector(t(d) %*% d))

      # Equation no. 14 from Smith & Abel, 1987.
      temp <- (1/2) *
        solve((((t(S) %*% P.orth_d) %*% W) %*% P.orth_d) %*% S) %*%
        ((((t(S) %*% P.orth_d) %*% W) %*% P.orth_d) %*% delta)

      R_s <- sqrt(sum(temp^2)) # Distance of origin from source
      error.temp <- delta - 2*(R_s*d) - 2*(S %*% temp)
      error <- sqrt(sum(error.temp^2)) # sum of squared error
      estimates[i,] <- data.frame(id = i, t(temp), error = error, eq = NA, stringsAsFactors = F)
    }else{
      #### Use normal form of localization (only works for 4 sensors)
      a <- 4 - 4*(t(d) %*% t(solve(S)) %*% solve(S) %*% d)
      b <- 2 * (t(d) %*% t(solve(S)) %*% solve(S) %*% delta) +
        2*(t(delta) %*% t(solve(S)) %*% solve(S) %*% d)
      C <- -(t(delta) %*% t(solve(S)) %*% solve(S) %*% delta)

      R_s.1 <- as.vector((-b + sqrt(b^2 - 4*a*C))/ (2*a))
      R_s.2 <- as.vector((-b - sqrt(b^2 - 4*a*C))/ (2*a))

      temp.1 <- 0.5 * solve(S) %*% (delta - (2*(R_s.1 * d)))
      temp.2 <- 0.5 * solve(S) %*% (delta - (2*(R_s.2 * d)))

      error.temp <- delta - 2*(R_s.1*d) - 2*(S %*% temp.1)
      error <- sqrt(sum(error.temp^2)) # calc sum of squared error

      estimates[((i-1)*2) + 1,] <- data.frame(id = i, t(temp.1), error = error, eq = '+', stringsAsFactors = F)
      estimates[((i-1)*2) + 2,] <- data.frame(id = i, t(temp.2), error = error, eq = '-', stringsAsFactors = F)
    }
  }

  ## Return hydrophones to real positions
  positions.real <- apply(estimates[,2:4],MARGIN = 1,FUN = function(x){
    x + hydrohpone.positions[1,]
  })
  positions.real <- do.call(what = 'rbind',args = positions.real)
  estimates[,2:4] <- positions.real
  return(estimates)
}

#' Find speed of sound from hydrophone array
#'
#' By minimizing error within the sheprical interpolation time-of-arrival-localization results and assuming the provided hydrophone positions are correct, this function will estimate the true speed of sound in your setup.
#' This is usefull for validating the speed of sound before processing with YAPS.
#'
#' @param c Your starting value for estimating the speed of sound.  Defaults to 1500m/s
#' @param sensorLocations A `n` x `3` array holding the xyz coordinates of the `n` hydrophones used in the setup.
#' @param detections A `source` x `hydrophone` array the detection times in seconds of the acoustic source.
#'
#' @return The estimated true speed of sound in m/s
#' @export
#'
findc <- function(c = 1500, sensorLocations, detections){
  message('Estimating speed of sound from TOA data...')

  ## Generate otimization function
  optFun <- function(c){
    loc <- TOA.localization(toa = detections,
                            hydrohpone.positions = sensorLocations,
                            c = c)
    if(nrow(sensorLocations) <= 4){
      idx.1 <- which(loc$eq == '-')
      idx.2 <- which(loc$eq == '+')
      if(sum(diff(loc$x[idx.1])^2 + diff(loc$y[idx.1])^2 + diff(loc$z[idx.1])^2) <
        sum(diff(loc$x[idx.2])^2 + diff(loc$y[idx.2])^2 + diff(loc$z[idx.2])^2)){
          return(sum(loc$error[idx.1]^2))
        }else{
          return(sum(loc$error[idx.2]^2))
        }
    }else{
      # Return sum of squared errors
      return(sum(loc$error^2))
    }
  }

  ## Run proccedure and return estimate
  vals <- optim(par = c, fn = optFun, method = 'Brent',
                lower = 1400, upper = 1600)
  return(vals)
}
