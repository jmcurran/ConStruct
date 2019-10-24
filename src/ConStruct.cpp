#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
double profileProb(IntegerVector::const_iterator prof, List freqs, int numLoci, double theta) {
  double p  = 1;

  for(int loc = 0; loc < numLoci; loc++){
    int i1 = 2 * loc;
    int a1 = prof[i1];
    int a2 = prof[i1 + 1];

    NumericVector f = as<NumericVector>(freqs[loc]);

    if(a1 == a2){
      double pA1 = f[a1 - 1];
      p *=  pA1 * (theta + (1 - theta) * pA1);
    }else{
      double pA1 = f[a1 - 1];
      double pA2 = f[a2 - 1];
      p *=  2 * pA1 * pA2 * (1 - theta);
    }
  }

  return p;
}

// [[Rcpp::export]]
double logLikelihood(double theta, IntegerVector Profiles, int numProfiles, int numLoci, List freqs) {
  double s = 0;
  IntegerVector::const_iterator it = Profiles.begin();

  for(int prof = 0; prof < numProfiles; prof++){
    s += log(profileProb(it, freqs, numLoci, theta));
    it += 2 * numLoci;
  }

  return s;
}

