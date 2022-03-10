

#include "KDA_intel.h"
#include "KDA_arm.h"

#include <Rcpp.h>
#include <cmath>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericVector calcKDA(const Rcpp::NumericVector& A)
{
  try {
    //convert abundances from A to Species
    size_t numspecies = A.size();
    std::vector< size_t > Abund(numspecies); //Abund = new int[numspecies];

    long double J = 0;
    for(size_t s = 0; s < numspecies; ++s) {
      if(s > Abund.size()) break;
      Abund[s] = A[s];
      J += Abund[s];
    }
    //call calcLogKDA
    std::vector<long double> K;


    if ( sizeof(long double) >= 16) {
      calcLogKDA(K, J, numspecies, Abund);
    } else {
      K = calcLogKDA_arm(numspecies, Abund);
    }

    int sizeofK = K.size();
    Rcpp::NumericVector out(sizeofK);
    for(int i = 0; i < sizeofK; ++i) {
      out[i] = K[i] + 4500.0 * logl(10);
    }

    return out;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
