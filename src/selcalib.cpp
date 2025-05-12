#include <TMB.hpp>

// Random Walk density ('huge' is optional)
template<class Type>
  Type RW_logdens(vector<Type> x, Type sd, int order=1){
    for(int i=0; i<order; i++) x = diff(x);
    return dnorm(x, Type(0), sd, true).sum();
  }
template<class Type>
  Type RW_logdens(vector<Type> x, Type sd, Type huge, int order=1){
    Type ans = 0;
    for(int i=0; i<order; i++) ans += dnorm(x(i), Type(0), huge, true);
    ans += RW_logdens(x, sd, order);
    return ans;
  }

// Main model
template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data and parameters */
  DATA_ARRAY(N);                // haul x size
  DATA_VECTOR(SweptArea);       // Sample area for each haul. length=nhaul
  DATA_FACTOR(group);           // length(group)=nhaul. Defines pair ids
  DATA_FACTOR(Gear);            // length(Gear)=nhaul. Defines the two gear types
  DATA_SCALAR(huge);            // "Infinite" standard deviation on log scale
  DATA_SCALAR(tiny);            // "Zero" standard deviation on log scale
  DATA_IVECTOR(rw_order);       // Order of logspectrum and loggear RW
  DATA_INTEGER(model_type);       // 0 = Poisson, 1 = ZIP, 2 = NB
  
  PARAMETER_ARRAY(logspectrum); // group x size. One spectrum for each pair id
  PARAMETER_ARRAY(nugget);      // haul x size. One nugget for each haul
  PARAMETER_ARRAY(residual);    // haul x size. One AR residual for each haul
  PARAMETER_VECTOR(loggear);    // 1 x size
  PARAMETER(logsd);             // sd of RW increments
  PARAMETER(phi);               // AR1(phi) contribution of residuals
  PARAMETER(logsdnug);          // sd of nugget effect
  PARAMETER(logsdres);          // sd of residuals
  PARAMETER(logsdGearRW);       // sd of increments in gear effect
  PARAMETER(logalpha);          // Exponent on the effect of SweptArea
  PARAMETER(logitpi);           // Only for ZIP
  PARAMETER(logtheta);          // Only for NB
  
  // Transpose for row-wise access
  array<Type> tnugget = nugget.transpose();
  array<Type> tresidual = residual.transpose();
  array<Type> tlogspectrum = logspectrum.transpose();
  array<Type> tN = N.transpose();
  
  int nhaul = N.dim[0];
  int nsize = N.dim[1];
  Type ans = 0;
  
  Type sd = exp(logsd);
  Type sdGearRW = exp(logsdGearRW);
  Type sdnug = exp(logsdnug);
  Type alpha = exp(logalpha);
  Type pi = Type(1) / (Type(1) + exp(-logitpi));
  Type theta_nb = exp(logtheta);
  
  // Random walk on log spectrum
  for (int i = 0; i < tlogspectrum.cols(); i++) {
    ans -= RW_logdens(vector<Type>(tlogspectrum.col(i)), sd, huge, rw_order(0));
  }
  
  // AR(1) residuals
  using namespace density;
  SCALE_t< AR1_t<N01<Type> > > nldens = SCALE(AR1(phi), exp(logsdres));
  for (int i = 0; i < nhaul; i++) {
    ans += nldens(tresidual.col(i));
  }
  
  // Nugget effect
  for (int j = 0; j < nsize; j++) {
    ans -= dnorm(vector<Type>(nugget.col(j)), Type(0), sdnug, true).sum();
  }
  
  // Gear effect RW
  ans -= RW_logdens(loggear, sdGearRW, huge, rw_order(1));
  
  // Likelihood
  for (int i = 0; i < nhaul; i++) {
    vector<Type> logintensity =
      tlogspectrum.col(group[i]) +
      tresidual.col(i) +
      alpha * log(SweptArea(i)) +
      tnugget.col(i);
    
    if (Gear(i) == 1) logintensity += loggear;
    else logintensity -= loggear;
    
    vector<Type> mu = exp(logintensity);
    vector<Type> y = tN.col(i);
    vector<Type> pois_ll = dpois(y, mu, true);
    
    switch(model_type) {

        case 1: // Poisson
          ans -= pois_ll.sum();
          break;

        case 2: // Zero-inflated Poisson
          for (int j = 0; j < nsize; j++) {
              if (y(j) == 0) {
                  ans -= log(pi + (1 - pi) * exp(pois_ll(j)));
              } else {
                  ans -= log(1 - pi) + pois_ll(j);
              }
          }
          break;

        case 3: // Negative Binomial
          for (int j = 0; j < nsize; j++) {
            Type prob = theta_nb / (theta_nb + mu(j));
            ans -= dnbinom(y(j), theta_nb, prob, true);
          }
          break;

    }
  }
  
  return ans;
}
