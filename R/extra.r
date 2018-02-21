#' Caclulate the circular mean
#'
#' @param rad A vector of radian values.
#' @param weight A vector of weights for radian values
#' @param verbose Can be usefull for debugging code.
#'
#' @return
#' @export
circular.mean <- function(rad,weight = c(), verbose = T){
  ## Get sample size
  n <- length(rad)

  if(n == 0){
    stop('Empty vector given')
  }

  ## Set weights to 1 if not specified
  if(missing(weight)){
    weight <- rep.int(1,times = n)
  }

  ## Set weights to 1 if not specified
  if(is.null(weight)){
    weight <- rep.int(1,times = n)
  }

  if(verbose == T)
    message('Mean weight: ',mean(weight))

  if(mean(weight) != 1){
    message('Average weighting is not equal to 1.  Adjusting weight values to correct this.')
    weight<- weight/mean(weight)
    message('Corrected mean weight: ',mean(weight))
  }

  ## Calculate the x and y components of the summed weighted vectors
  sum.y <- sum(weight*sin(rad))
  sum.x <- sum(weight*cos(rad))

  ## Get the mean resultant length
  rbar <- (sqrt(sum.x^2 + sum.y^2))/sum(weight)

  theta <- atan2(sum.y,sum.x)

  return(list(theta = theta,rbar = rbar))
}
