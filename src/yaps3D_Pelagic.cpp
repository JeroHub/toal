#include <TMB.hpp>
using namespace density;
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Model input (data)
  DATA_ARRAY(H);			// Position of hydros
  DATA_ARRAY(toa);   		// Time of arrival at hydro. One row per buoy, one column per ping
  DATA_INTEGER(nh); // number of hydrohpones
  DATA_INTEGER(np); // number of pings

  /***************************
  * Parameters
  ***************************/
  // Random
  //--------------------------
  PARAMETER_ARRAY(XYZ);	    // Estimated random position at time of ping
  PARAMETER_VECTOR(top);		// Estimated random times of pings
  PARAMETER_ARRAY(dl);      // Residual latency for each detection

  // Fixed
  //--------------------------
  DATA_INTEGER(c); // Speed of sound
  PARAMETER(logSigma_bi);
  Type sigma_bi = exp(logSigma_bi);

  // Hydrophone latency distribution
  PARAMETER(logSigma_dl);

  // TOA distribution
  PARAMETER(logSigma_toa); // Time of arrival SD (for each hydrophone)
  // PARAMETER(logScale_toa);		// scale-parameter for t-dist
  // PARAMETER(log_t_part);		// t-part of mixture model

  PARAMETER(logSigma_xyz);

  /*****************************
  * Helper Variables
  *
  * Variables for temporary storage to make
  * life easier.
  * Try to keep to a minimum
  *****************************/
  array<Type> mu_toa(np,nh);  // Estimated toa
  array<Type> dist(np,nh);	 // dist-matrix

  // SD for hydrophone latency errors
  Type sigma_dl = exp(logSigma_dl);
  // vector<Type> sigma_dl(nh);
  // for(int i = 0; i < nh; ++i){
  //   sigma_dl(i) = exp(logSigma_dl(i));
  // }

  // Mixture of guassian and t-distriubtions for TOA
  Type sigma_toa = exp(logSigma_toa);
  // Type scale_toa = exp(logScale_toa);
  // Type t_part = exp(log_t_part);
  // Type G_part = Type(1.0) - t_part;

  // Displacement
  Type sigma_xyz = exp(logSigma_xyz);
  Type velocity;

  /***************************************************
  * Start run
  ***************************************************/
  // Negative log liklihood
  // Use the non-parallel intialization for debugging
  // Type nll = 0;
  parallel_accumulator<Type> nll(this);

  /*************************************************
  * Time of arrival
  * Assumes a mixture of Guassian and t error distributions
  * for residuals of ping detections (eps)
  *************************************************/
  for(int i=0; i<np; ++i){  //iterate pings
    for(int h=0; h<nh; ++h){ //iterate hydros
      if(toa(i,h) != -9999){ //ignore NA's...
        // Calculate Distance from hydrophone
        dist(i,h) = sqrt(
            (H(h,0)-XYZ(i,0))*(H(h,0)-XYZ(i,0)) +
            (H(h,1)-XYZ(i,1))*(H(h,1)-XYZ(i,1)) +
            (H(h,2)-XYZ(i,2))*(H(h,2)-XYZ(i,2)));

        // Expected time of arrival
        mu_toa(i,h) = top(i) + (dist(i,h)/c) + dl(i,h);

        // Residual time of arrival (real - expected)
        Type eps = toa(i,h) - mu_toa(i,h);

        // Negative log liklihood from time of arrival
        // This is a mixture of gaussian and scaled t-distribution.
        // Does the t-distribution help with smaller sample sizes?
        // nll -= log( G_part * dnorm(eps, Type(0), sigma_toa, false) + //Gaussian part
        //   t_part * dt(eps/scale_toa, Type(3.0), false) / scale_toa );	//t part

        // Distribution of residuals
        nll -= dnorm(eps, Type(0), sigma_toa, true);
      }
    }
  }

  /*************************************************
  * Spatial location
  * Assumes a gamma error distribution of swim displacement
  * along each axis between pings.
  *************************************************/
  for(int i=1; i<np; ++i){
    // displacement liklihood
    nll -= dnorm(XYZ(i,0), XYZ(i-1,0),
                 sqrt(sigma_xyz*(top(i) - top(i-1))), true);
    nll -= dnorm(XYZ(i,1), XYZ(i-1,1),
                 sqrt(sigma_xyz*(top(i) - top(i-1))), true);
    nll -= dnorm(XYZ(i,2), XYZ(i-1,2),
                 sqrt(sigma_xyz*(top(i) - top(i-1))), true);


    // Take sqrt of dist to transfrom into a normal distribution
    // The mean should be linearly related to the SD if x, y, and z displacement
    // is normally distributed and centered around 0
    // velocity = sqrt(
    //   (XYZ(i,0) - XYZ(i-1,0))*(XYZ(i,0) - XYZ(i-1,0)) +
    //     (XYZ(i,1) - XYZ(i-1,1))*(XYZ(i,1) - XYZ(i-1,1)) +
    //     (XYZ(i,2) - XYZ(i-1,2))*(XYZ(i,2) - XYZ(i-1,2)))/
    //       (top(i) - top(i-1));
    //
    // nll -= dnorm(sqrt(velocity), sigma_xyz*4.508148, sigma_xyz, true);
  }

  /*************************************************
   * Hydrophone detection latency (+/-)
   * Caused by: sample rate limitations and
   * hydrophone movement,
   *************************************************/

  for(int i = 0; i<np; i++){
    for(int h = 0; h<nh; h++){
      nll -= dnorm(dl(i,h), Type(0), sigma_dl, true);
    }
  }

  /*************************************************
   * Inter-ping-interval restrictions
   * Assumes a constant ping rate, with slight fluctuations
   *************************************************/

  //burst interval component
  for(int i = 2; i < np; ++i)
  {
      nll -= dnorm(top(i)-2*top(i-1)+top(i-2), Type(0),sigma_bi, true);
  }

  return nll;
}
