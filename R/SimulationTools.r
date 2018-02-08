#' Convert xyz positions to time of arrival values for a hydrophone array setup
#'
#' This can be used to either calculate time-of-arrival information from a given set of tag and hydrophone positions.
#' The parameters `reciever.fs`, `hydrophoneMovement.sd`, and `p.detect` can be used to simulated sources of error from the detection systems sample rate, hydrophone movement, and probabolity of detection, respectively.
#'
#' The parameter `hydrophoneMovement.sd` will simulate spatial error by applying a guassian normal error value to the distance between the hydrophone and the tag.
#' This is intended to simulate random hydrophone movement in three dimentions which way be the result of wave action (in particular for hydrophones mounted on floating structures).
#'
#' @param hydrophone.positions An `n`-by-3 matrix or data.frame, holding the x, y, and z positions, in meters, of the hydrophones where `n`` is the number of hydrophones in the array.
#' @param tag.positions An `n`-by-4 matrix or dataframe holding the x, y, z positions and time in seconds of the tag pulses where `n` is the number of pulses.
#' @param c Speed of sound through the medium (meters per second).
#'
#' @return
#' @export
xyz2toa <- function(hydrophone.positions,
                    tag.positions,
                    c = 1500){
  ## Get sample size
  n.hydro <- nrow(hydrophone.positions) # Number of sensors
  n.pings <- nrow(tag.positions)

  ## Helper function for calculating distance between points
  distance <- function(a,b){
    ## Calc distance between two locations
    dif <- abs(a - b)
    return(sqrt(sum(dif^2)))
  }

  #############################
  ## Initalize toa array
  #############################
  toa <- matrix(ncol = n.hydro,nrow = n.pings)
  # Calculate real toa for each ping and reciever
  for(i in 1:n.hydro){
    browser()
    toa[,i] <- apply(
      X = tag.positions,
      MARGIN = 1,
      FUN = function(x){
        browser()
        return(distance(x[1:3],hydrophone.positions[i,1:3]) + x[4])
      })
  }

  ## Label and return matrix
  rownames(toa) <- paste0('Ping ',1:nrow(tag.positions))
  colnames(toa) <- paste0('Hydrophone ',1:N)
  return(toa)
}

#' Simulate error on an array of time-of-arrival values
#'
#' The parameter `hydrophoneMovement.sd` will simulate spatial error by applying a guassian normal error value to the distance between the hydrophone and the tag.
#' This is intended to simulate random hydrophone movement in three dimentions which way be the result of wave action (in particular for hydrophones mounted on floating structures).
#'
#'
#' @param toa An `n`-by-`m` matrix or data.frame, holding the time-of-arrival values, in seconds.
#' The number of pings is `n`  while `m` is the number of hydrophones in the array.
#' @param tag.positions An `n`-by-4 matrix or dataframe holding the x, y, z positions and time in seconds of the tag pulses where `n` is the number of pulses.
#' @param c Speed of sound through the medium (meters per second).
#' @param reciever.fs The sample rate of the reciever (Hz) used to simulate toa error.
#' @param hydrophoneMovement.sd The standard deviation of the 0 centered Guassian distribution used to simulate toa error due to hydrophone movement.
#' A single value to be applied to all hydrophones, or a vector which equals the length of the number of hydrophones.
#' @param p.detect The probability (between 0 and 1) that a hydrophone will detect a pulse.
#' A single value to apply to all hydrophones, or a vector of values the same length as the number of hydrophones.
#'
#' @return
#' @export
toaSimError <- function(toa, reciever.fs = c(),
                        hydrophoneMovement.sd = c(), p.detect = c()){

  #############################################################
  ## Error due to hydrophone movement (gaussian distributed)
  #############################################################
  if(length(hydrophoneMovement.sd) > 0){
    if(length(hydrophoneMovement.sd) == 1){
      ## Single error value for all hydrophones
      error.hydroMovement <- rnorm(n = length(toa),
                                   mean = 0, sd = hydrophoneMovement.sd)
      toa <- toa + (error.hydroMovement/c)
    }else if(length(hydrophoneMovement.sd) == n.hydro){
      # Unique error values for each hydrophone
      for(i in 1:n.hydro){
        error.hydroMovement <- rnorm(n = length(toa[,i]),
                                     mean = 0, sd = hydrophoneMovement.sd)
        toa[,i] <- toa[,i] + (error.hydroMovement/c)
      }
    }else{
      message('Parameter hydrophoneMovement.sd must be a length equal to 1 or the number of hydrophones in the reciever array (', n.hydro,').')
    }
  }

  ###########################################
  ## Apply error due to minimum sample rate
  ###########################################
  if(length(reciever.fs) > 0){
    toa <- round(toa*reciever.fs,digits = 0)*(1/reciever.fs)
  }

  #############################
  ## Apply missed detections
  #############################
  if(length(p.detect) > 0){
    if(length(p.detect) == 1){
      ## Single error value for all hydrophones
      ## randomly assign NA values to each detection
      error.detect <- rbinom(n = length(toa),
                             prob = p.detect,
                             size = 1)
      toa[which(error.detect)] <- NA
    }else if(length(p.detect) == n.hydro){
      # Unique error values for each hydrophone
      for(i in 1:n.hydro){
        error.detect <- rbinom(n = length(toa[,i]),
                               prob = p.detect,
                               size = 1)
        toa[which(error.detect),i] <- NA
      }
    }else{
      message('Parameter p.detect must be a length equal to 1 or the number of hydrophones in the reciever array (', n.hydro,').')
    }
  }
  return(toa)
}
