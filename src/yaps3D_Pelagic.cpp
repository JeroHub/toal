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
	* Parameters (Fixed)
	***************************/
	// Speed of sound
	PARAMETER(c);

	// Time of arrival
	PARAMETER(logSigma_toa);
	Type sigma_toa = exp(logSigma_toa);

	/***************************
	* Parameters (Random)
	* These must be Guassian distributed
	*  (for Laplace approximation)
	***************************/
	// Spatial location at time of ping
	PARAMETER_ARRAY(XYZ);	    // Estimated random position at time of ping
	PARAMETER(logD_xy);    		// Log SD of XY movement/unit time
	Type D_xy = exp(logD_xy);
  PARAMETER(logD_z);        // Log SD of Z movement/unit time
	Type D_z = exp(logD_z);

	// Time of pings
	PARAMETER_VECTOR(top);		// Estimated random time of pings
	PARAMETER(logSigma_bi);		// Sigma for burst interval
	Type sigma_bi = exp(logSigma_bi);

	// Hydrophone movement (along axis facing tag)
	DATA_ARRAY(tag_movement);
	PARAMETER(logD_xy);    		// Log SD of XY movement/unit time
	Type D_xy = exp(logD_xy);

	/********************************
	* Distribution parameters
	* For mixing of distributions
	********************************/
	// Mixed Guassian/t for toa estimates
	PARAMETER(logScale_toa);		// scale-parameter for t-dist
	Type scale = exp(logScale_toa);
	PARAMETER(log_t_part);		// t-part of mixture model
	Type t_part = exp(log_t_part);
	Type G_part = Type(1.0) - t_part; // Gaussian part of mixture model

	/*****************************
	* Helper Variables
	*****************************/
	array<Type> mu_toa(nh,np);  // Estimated toa
	array<Type> dist(nh,np);	 // dist-matrix

	// Type nll (negative log liklihood);
	// Lower values is better fitting model
	parallel_accumulator<Type> nll(this);

/*************************************************
 * Time of arrival
 * Assumes a mixture of Guassian and t error distributions
 * for residuals of ping detections (eps)
 *************************************************/
	for(int i=0; i<np; ++i) //iterate pings
	{
		for(int h=0; h<nh; ++h){ //iterate hydros
			if(toa(h,i) != -9999){ //ignore NA's...
			  // Calculate Distance from hydrohone
				dist(h,i) = sqrt(
				    (H(h,0)-XYZ(i,0))*(H(h,0)-XYZ(i,0)) +
				    (H(h,1)-XYZ(i,1))*(H(h,1)-XYZ(i,1)) +
				    (H(h,2)-XYZ(i,2))*(H(h,2)-XYZ(i,2))
			  );

			  // Calc expected toa from distance and estimated ping time
				mu_toa(h,i) = top(i) + (dist(h,i)/c);

				// Residual time of arrival (real - expected)
				Type eps = toa(h,i) - mu_toa(h,i);

				// Negative log liklihood from time of arrival
				// This is a mixture of gaussian and scaled t-distribution.
				// Does the t-distribution help with smaller sample sizes?
				nll -= log( G_part * dnorm(eps, Type(0),sigma_toa,false) + //Gaussian part
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
	for(int i=0; i<np; ++i)
	{
		if(i == 0) {
		 	nll -= dnorm(XYZ(0,0),Type(0),Type(1000),true);
			nll -= dnorm(XYZ(0,1),Type(0),Type(1000),true);
			nll -= dnorm(XYZ(0,2),Type(0),Type(1000),true);
		} else {
			nll -= dnorm(XYZ(i,0),XYZ(i-1,0),
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
 * Pulse intervals
 *************************************************/
	// Burst interval component
	// Burst interavals (time between transmission emissions) are random guassian
	// This allows for variations in tag intervals due to temperature
	// Also compensates for sample rate limitations?
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
