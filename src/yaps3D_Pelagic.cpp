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
  PARAMETER(c); // Speed of sound
  PARAMETER(logSigma_bi);
  Type sigma_bi = exp(logSigma_bi);

  // Hydrophone latency distribution
  PARAMETER_VECTOR(logSigma_dl);

  // XYZ distribution
  PARAMETER(logD_xy);    		// Log SD of XY movement/unit time
  PARAMETER(logD_z);        // Log SD of Z movement/unit time

  // TOA distribution
  PARAMETER(logSigma_toa); // Time of arrival SD (for each hydrophone)
  PARAMETER(logScale_toa);		// scale-parameter for t-dist
  PARAMETER(log_t_part);		// t-part of mixture model

  /*****************************
  * Helper Variables
  *****************************/
  array<Type> mu_toa(np,nh);  // Estimated toa
  array<Type> dist(np,nh);	 // dist-matrix

  vector<Type> sigma_dl(nh);// SD for hydrophone latency errors
  for(int i = 0; i < nh; ++i){
    sigma_dl(i) = exp(logSigma_dl(i));
  }

  Type sigma_toa = exp(logSigma_toa);
  Type D_xy = exp(logD_xy);
  Type D_z = exp(logD_z);
  Type scale_toa = exp(logScale_toa);
  Type t_part = exp(log_t_part);
  Type G_part = Type(1.0) - t_part; // Gaussian part of mixture model
  Type x1, y1, z1, mag, x2, y2, z2;  // For calculating delta movement vector

  // Type nll (negative log liklihood);
  // Lower values is better fitting model
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
        nll -= log( G_part * dnorm(eps, Type(0), sigma_toa, false) + //Gaussian part
          t_part * dt(eps/scale_toa, Type(3.0), false) / scale_toa );	//t part
      }
    }
  }

  /*************************************************
  * Spatial location
  * Assumes a guassian error distribution of distances
  * traveled along x, y, and z axes between pings.
  * SD is scaled by inter pulse interval
  *************************************************/
  for(int i=0; i<np; ++i){
    if(i == 0) {
      nll -= dnorm(XYZ(0,0),Type(0),Type(1000),true);
      nll -= dnorm(XYZ(0,1),Type(0),Type(1000),true);
      nll -= dnorm(XYZ(0,2),Type(0),Type(1000),true);
    } else {
      nll -= dnorm(XYZ(i,0),
                   XYZ(i-1,0),
                   sqrt(2*D_xy*(top(i) - top(i-1))),
                   true);
      nll -= dnorm(XYZ(i,1),
                   XYZ(i-1,1),
                   sqrt(2*D_xy*(top(i) - top(i-1))),
                   true);
      nll -= dnorm(XYZ(i,2),
                   XYZ(i-1,2),
                   sqrt(2*D_z*(top(i) - top(i-1))),
                   true);
    }
  }

  /*************************************************
   * Hydrophone detection latency (+/-)
   * Caused by: sample rate limitations and
   * hydrophone movement,
   *************************************************/

  for(int i = 0; i<np; i++){
    for(int h = 0; h<nh; h++){
      nll -= dnorm(dl(i,h), Type(0), sigma_dl(h), true);
    }
  }

  /*************************************************
   * Inter-ping-interval restrictions
   * Assumes a constant ping rate, with slight fluctuations
   *************************************************/

  //burst interval component
  for(int i = 2; i < np; ++i)
  {
    if(i == 0) {
      nll -= dnorm(top(0),Type(0.0),Type(4.0),true);
    } else if (i == 1){
      nll -= dnorm(top(1),Type(2.0),Type(4.0),true);
    } else {
      nll -= dnorm(top(i)-2*top(i-1)+top(i-2), Type(0),sigma_bi, true);
    }
  }

  return nll;
}
