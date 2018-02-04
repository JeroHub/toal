#include <TMB.hpp>
using namespace density;
template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_ARRAY(H);			// Position of hydros
	DATA_ARRAY(toa);   		// Time of arrival at hydro. One row per buoy, one column per ping
	DATA_INTEGER(nh); // number of hydrohpones
	DATA_INTEGER(np); // number of pings

	PARAMETER_ARRAY(XY);	//Position at time of ping
	PARAMETER_VECTOR(top);		// Estimated time of pings

	PARAMETER(logD_xy);    		// Diffusivity of fish movement distances (i.e. log variance)
	Type D_xy = exp(logD_xy);
	PARAMETER(logD_z);
	Type D_xy = exp(logD_z);

	PARAMETER(logSigma_bi);		// Sigma for burst interval
	Type sigma_bi = exp(logSigma_bi);

	PARAMETER(ss);    		// Speed of sound

	PARAMETER(logSigma_toa);	// Sigma TimeOfArrival
	Type sigma_toa = exp(logSigma_toa);

	PARAMETER(logScale);		// scale-parameter for t-dist
	Type scale = exp(logScale);

	PARAMETER(log_t_part);		// t-part of mixture model
	Type t_part = exp(log_t_part);
	Type G_part = Type(1.0) - t_part; // Gaussian part of mixture model

	array<Type> mu_toa(nh,np);  // mu-matrix
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
				    (H(h,0)-XY(i,0))*(H(h,0)-XY(i,0)) +
				    (H(h,1)-XY(i,1))*(H(h,1)-XY(i,1)) +
				    (H(h,2)-XY(i,2))*(H(h,2)-XY(i,2))
			  );

			  // Calc expected toa from distance and estimated ping time
				mu_toa(h,i) = top(i) + (dist(h,i)/ss);

				// Residual time of arrival (real - expected)
				Type eps = toa(h,i) - mu_toa(h,i);

				// Negative log liklihood from time of arrival
				// This is a mixture of gaussian and scaled t-distribution.
				// Does the t-distribution help with smaller sample sizes?
				nll -= log( G_part * dnorm(eps, Type(0),sigma_toa,false) + //Gaussian part
							t_part * dt(eps/scale, Type(3.0), false) / scale );	//t part
			}
		}
	}

/*************************************************
 * Spatial location
 * Assumes a guassian error distribution of distances
 * traveled along x, y, and z axes between pings
 *************************************************/
	for(int i=0; i<np; ++i)
	{
		if(i == 0) {
		 	nll -= dnorm(XY(0,0),Type(0),Type(1000),true);
			nll -= dnorm(XY(0,1),Type(0),Type(1000),true);
			nll -= dnorm(XY(0,2),Type(0),Type(1000),true);
		} else {
			nll -= dnorm(XY(i,0),
                XY(i-1,0),
                sqrt(2*D_xy*(top(i) - top(i-1))),
                true);
			nll -= dnorm(XY(i,1),
                XY(i-1,1),
                sqrt(2*D_xy*(top(i) - top(i-1))),
                true);
			nll -= dnorm(XY(i,2),
                XY(i-1,2),
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
