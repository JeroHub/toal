#include <TMB.hpp>
using namespace density;
template<class Type>
Type objective_function<Type>::operator() ()
{

  // -999: Code for missing data

  /**** Input ****/
  DATA_IARRAY(ID);			// [Trial, Group, Member] array for individual IDs
  DATA_ARRAY(OxygenGroup);// [Trial, Group, Speed, Measure] array for group O2 measures
  DATA_ARRAY(Temperature);// [Trial, Group, Speed, Measure] array for group Temp measures

  DATA_VECTOR(Speed);// Speeds in m/s
  DATA_VECTOR(Length);// Speeds in m/s

  DATA_INTEGER(nTrials);
  DATA_INTEGER(nGroups);
  DATA_INTEGER(nMembers);

  DATA_INTEGER(nSpeeds);
  DATA_INTEGER(nMeasures);

  /**** Parameters ****/
  // ---- Random ----
  PARAMETER_VECTOR(BaseMetabBias); // [ID] Individual bias (random effect) for O2 consumption

  // ---- Fixed ----
  // Covariates
  PARAMETER(BaseMetab); // base metabolism: O2 consumption @ resting per length (g)
  PARAMETER(LengthMetab); // base metabolism: O2 consumption @ resting per length (g)
  PARAMETER(SwimMetab); // O2 consumption increase per speed (m/s)
  PARAMETER(TempEffect); // O2 consumption multiplier as result of temerature

  // Variance
  PARAMETER(logSigma_Residuals); // Log SD for O2 consumption per weight (g)
  PARAMETER(logSigma_BaseMetabBias); // Log SD for O2 consumption per weight (g)

  /**** Helper Variables ****/
  Type Sigma_Residuals = exp(logSigma_Residuals);
  Type Sigma_BaseMetabBias = exp(logSigma_BaseMetabBias);

  Type groupO2;
  int id;

  /**** StaRt Run ****/
  Type nll = 0;
  //parallel_accumulator<Type> nll(this);


  /*************************************************
  * Group O2
  * Assumes a mixture of Guassian and t error distributions
  * for residuals of ping detections (eps)
  *************************************************/
  for(int trial=0; trial<nTrials; ++trial){  //Iterate trials
    for(int group=0; group<nGroups; ++group){  //Iterate groups
      for(int speed=0; speed<nSpeeds; ++speed){  //Iterate trials
        for(int measure=0; measure<nMeasures; ++measure){  //Iterate measures (i.e. time intervals)

          // Identify individuals in group

          // Calculate predicted group O2 consumption
          groupO2 = 0;
          for(int member=0; member < nMembers; ++member){
            // check for NA data (missing group members or O2 measures
            if (ID(trial,group,member) != -999 &&
                OxygenGroup(trial,group,speed,measure) != -999 &&
                Temperature(trial,group,speed,measure) != -999){

              // O2 = base metabolism + length + individual random effect + speed effect
              groupO2 = BaseMetab + LengthMetab*Length(id) + BaseMetabBias(id) + (SwimMetab*Speed(speed));

              // Add temperature effects (linear relationship between temperature and O2)
              groupO2 = groupO2 * (TempEffect*Temperature(trial,group,speed,measure));

              // Calculate individual bias liklihood (Random effects)
              nll -= dnorm(BaseMetabBias(id), Type(0), Sigma_BaseMetabBias, true);
            }
          }

          // Calculate residual liklihood (Fixed effects)
          nll -= dnorm(groupO2 - OxygenGroup(trial,group,speed,measure), Type(0),
                       Sigma_Residuals, true);
        }
      }
    }
  }

  return nll;
}
