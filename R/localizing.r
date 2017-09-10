#' Time-of-arrival difference localization
#'
#' Localize sound sources using an array of 4 or more hydrophones.
#' Four hydrophones will result in 1 degree of ambiguity and two position estimates, while more hydrophones will result in a single estimate.
#' This is an implementation by the spherical interpolation methods described by Schau & Robinson (1987) and Smith & Abel (1987) which provides a closed form solution for estimating source locations from the relative differences between time-of-arrival at each hydrophone.
#'
#' The parameter \code{toa} must be a 4 column matrix or data.frame where each column represents a unique hydrophone in the array and each row represents a detection.
#' The values in this parameter represent the time of arrival in seconds at each hydrophone.
#' The parameter \code{hydrophone.positions} should be a N by 3 matrix or data.frame, where each row represents a hydrophone (in the same order as appearing in the columns of \code{toa}) and the three columns are the x, y, and z positions in meters.
#'
#' @param toa data.frame or matrix holding times of arrival for each hydrophone  See details for formatting requirements.
#' @param hydrophone.positions data.frame or matrix holding information about hydrophone positions.  See details for formatting.
#' @param c Speed of sound in the medium.  The default value represents a general approximation of sea water.
#'
#' @return
#' @export
#'
#' @examples
#'
#' @references Smith, J. O., Abel, J. S. 1987. Closed-Form Least-Squares Source Location Estimation from Range-Difference Measurements.  IEEE TRANSACTIONS ON ACOUSTICS, SPEECH, AND SIGNAL PROCESSING. 35(12):1661-1669.
#'
#' Schau, H. C., Robinson, A. Z. 1987. Passive Source Localization Employing Intersecting Spherical Surfaces from Time-of-Arrival Differences.  IEEE TRANSACTIONS ON ACOUSTICS, SPEECH, AND SIGNAL PROCESSING. 35(8):1223-1225
TOA_localization <- function(toa, hydrohpone.positions, c = 1500){

  ## Nested function for calculating distance between points in 3d space
  dist <- function(a,b){
    ## Calc distance between two locations
    dif <- abs(a - b)
    return(sqrt(sum(dif^2)))
  }

  ## Define S: sensor locations, use first sensor as origin
  # x,y,z position of sensors 2-4 (sensor = receivers)
  S  = as.matrix(sensor.positions[2:N,c('x','y','z')])

  ## Define d: Range distance differences distances between sensors and reference sensor
  ## Crete matrix of distances between sensors i,j
  d <- matrix(nrow = N-1,ncol = 1)
  for(i in 1:(N-1)){
    TOA.diff <- sensorLatency[i+1] - sensorLatency[1]
    range.diff <- TOA.diff*c # distance differences
    d[i,1] <- range.diff
  }
  ## Add row names to matrix, defining sensors i and j
  row.names(d) <- paste0(as.character(2:N),
                         ',',
                         as.character(rep.int(1,times = N-1)))
  colnames(d) <- 'RD between i,j'

  ## Define R: Distances between Reference sensor (id: 1) and other sensors
  R <- matrix(nrow = N-1,ncol = 1)
  for(i in 1:(N-1)){
    R[i,1] <- sqrt(sum(S[i,]^2))
  }
  ## Add row names to matrix, defining sensors i and j
  row.names(R) <- as.character(2:N)
  colnames(R) <- 'Distance from sensor i to reference sensor (id: 1)'

  ## delta: R - d
  delta <- R^2 - d^2
  ## 'Weighting matrix.
  # Set all values to 1 and assume all sesnors are equally reliable
  W = diag(x = rep.int(1,times = N-1))

  ## Identity matrix
  I = diag(x = rep.int(1,times = N-1))

  ## Orthogonal compliment of P, with respect to d
  P.orth_d <- I - ((d %*% t(d)) / as.vector(t(d) %*% d))

  # Equation no. 14 from Smith & Abel, 1987.
  x_source = (1/2) *
    ginv((((t(S) %*% P.orth_d) %*% W) %*% P.orth_d) %*% S) %*%
    ((((t(S) %*% P.orth_d) %*% W) %*% P.orth_d) %*% delta)

  return(x_source)
}
