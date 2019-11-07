#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
double profileProb(IntegerVector::const_iterator prof, List freqs, int numLoci, double theta, double r = 0) {
  double p  = 1;

  if(r != 0){
    for(int loc = 0; loc < numLoci; loc++){
      int i1 = 2 * loc;
      int a1 = prof[i1];
      int a2 = prof[i1 + 1];

      NumericVector f = as<NumericVector>(freqs[loc]);

      if(a1 == a2){
        double pA1 = f[a1 - 1];
        p *=  pA1 * (r + (1 - r)* (theta + (1 - theta) * pA1));
      }else{
        double pA1 = f[a1 - 1];
        double pA2 = f[a2 - 1];
        p *=  2 * pA1 * pA2 * (1 - theta) * (1 - r);
      }
    }
  }else{
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

// [[Rcpp::export]]
double logLikelihoodCosang(double theta, double C, double r, IntegerVector Profiles, int numProfiles, int numLoci, List freqs) {
  double s = 0;
  IntegerVector::const_iterator it = Profiles.begin();

  for(int prof = 0; prof < numProfiles; prof++){
    s += log((1 - C) * profileProb(it, freqs, numLoci, theta) + C * profileProb(it, freqs, numLoci, theta, r)) ;
    it += 2 * numLoci;
  }

  return s;
}

// [[Rcpp::export]]
NumericMatrix logLikelihoodCosangMat(NumericVector theta, NumericVector C, double r,
                                  IntegerVector Profiles, int numProfiles, int numLoci, List freqs) {
  NumericMatrix result(theta.size(), C.size());
  IntegerVector::const_iterator it;
  double F,pi;
  double j = 0;
  double s = 0;

  NumericMatrix probSubstruct(theta.size(), numProfiles);
  NumericMatrix probCosang(theta.size(), numProfiles);

  for(int row = 0; row < theta.size(); row++){
    F = theta[row];
    it = Profiles.begin();

    for(int prof = 0; prof < numProfiles; prof++){
      probSubstruct(row, prof) = profileProb(it, freqs, numLoci, F);
      probCosang(row, prof) = profileProb(it, freqs, numLoci, F, r);
      it += 2 * numLoci;
    }
  }

  for(int row = 0; row < theta.size(); row++){
    for(int col = 0; col < C.size(); col++){
      pi = C[col];
      s = 0;

      for(int prof = 0; prof < numProfiles; prof++){
        s += log((1 - pi) * probSubstruct(row, prof) + pi * probCosang(row, prof));
      }

      result(row, col) = s;
    }
  }

  return NumericMatrix(theta.size(), C.size(), result.begin());
}



