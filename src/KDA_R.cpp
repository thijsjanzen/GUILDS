// Copyright 2015 - 2025 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include "KDA_arm.h"          // NOLINT [build/include_subdir]

#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>   // qsort
#include <limits.h>
#include <math.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>



using namespace std;         // NOLINT [build/namespaces]
using namespace Rcpp;        // NOLINT [build/namespaces]
//// Base of this code by
//   J. Chave and F. Jabot
///  last update 05-23-2008

std::vector<long double> calcLogKDA(int numspecies,
                                    std::vector<int> Abund)  {
  std::vector<long double> K;

  if (Abund.size() < 1) return K;

  sort(Abund.begin(), Abund.end());


  long double J = 0;
  long double SPP = numspecies;
  for (int s = 0; s < SPP; s++) {
    J += Abund[s];
  }

  int MaxA = 0;
  MaxA = Abund[ Abund.size() - 1];    // MaxA = Abund[(int)SPP-1];

  // abundance distribution
  // int *Phi = new int[MaxA+1]; for(s=0;s<=MaxA;s++) Phi[s]=0;
  std::vector<int> Phi(MaxA + 1, 0);

  for (int s = 0; s < SPP; s++) {
    Phi[Abund[s]]++;
  }

  // Number of distinct abundances
  int NDA = 0;
  for (int s = 0; s <= MaxA; s++) {
    if (Phi[s] > 0) {
      NDA++;
    }
  }

  // cerr << "Start computing Stirling numbers ...\n";
  // FIRST STAGE: compute the Stirling numbers
  // I use the relationship S(n,m)=S(n-1,m-1)-(n-1)S(n-1,m)
  // where n >= m >= 1
  // The relation of interest is sum_m S(n,m)/S(n,1)*S(m,1) * x^m
  // defined to be equal to sum_m T(n,m) * x^m
  // The recurrence relation on T(n,m) is
  // T(n,m)= T(n-1,m) + T(n-1,m-1)*(m-1)/(n-1)


  std::vector<int> f(NDA, 0);   // int *f = new int[NDA];
  std::vector<int> g(NDA, 0);   // int *g = new int[NDA];
  int i = 0;
  // for(s=0;s<NDA;s++) {f[s]=0;g[s]=0;}
  for (int n = 0; n <= MaxA; n++) {
    if (Phi[n] > 0) {
      f[i] = Phi[n];
      g[i] = n;
      i++;
    }
  }
  // long double **T= new long double*[NDA];
  // T(n,m) just for the n which are useful
  // T[0] = new long double[g[0]+1];


  std::vector< long double > fill;
  std::vector< std::vector< long double > > T(NDA, fill);

  std::vector<long double> T0(g[0] + 1);
  T[0] = T0;
  T[0][0] = 0; T[0][1] = 1;

  if (g[0] != 1) {
    //  long double *lS2 = new long double[g[0]+1];
    std::vector<long double> lS2(g[0] + 1);
    lS2[0] = 0; lS2[1] = 1;
    for (int n = 2; n <= g[0]; n++) {
      // long double *lS1 = new long double[n+1];
      std::vector<long double> lS1(n + 1);
      for (int im = 0; im <= n - 1; im++) {
        lS1[im] = lS2[im];
      }
      lS1[n] = 0;
      for (int im = 2; im <= n; im++) {
        lS2[im] = lS1[im] + lS1[im - 1] * (im - 1) / (n - 1);
      }
    }
    for (int im = 2; im <= g[0]; im++) {
      T[0][im] = lS2[im];
    }
  }


  for (int in = 1; in < i; in++) {
    // T[in]= new long double[g[in]+1];
    std::vector<long double> Tin(g[in] + 1);
    T[in] = Tin;
    T[in][0] = 0; T[in][1] = 1;

    // long double *lS2 = new long double[g[in]+1];
    std::vector<long double> lS2(g[in] + 1);
    for (int im = 0; im <= g[in - 1]; im++) {
      lS2[im] = T[in - 1][im];
    }
    for (int n = g[in - 1] + 1; n <= g[in]; n++) {
      // long double *lS1 = new long double[n+1];
      std::vector< long double > lS1(n + 1);
      for (int im = 0; im <= n - 1; im++) {
        lS1[im] = lS2[im];
      }
      lS1[n] = 0;
      for (int im = 2; im <= n; im++) {
        lS2[im] = lS1[im] + lS1[im - 1] * (im - 1) / (n - 1);
      }
    }
    for (int im = 2; im <= g[in]; im++) {
      T[in][im] = lS2[im];
    }
  }
  // After this stage we have stored in T[i][m] T(g[i],m)
  // with T(n,m) = S(n,m)*S(m,1)/S(n,1) for i>0

  // cerr << "Start computing ln(K(D,A)) ...\n";
  // SECOND STAGE: compute the K(D,A)
  // I follow Etienne's route. Compute the product of polynomials
  // of length J
  K.clear();
  K.resize(J + 2, 0.0);   // K = new long double[int(J)+1];

  // long double *poly2 = new long double[int(J)+1];
  std::vector<long double> poly2(J + 1, 0.0);
  K[0] = 1;
  int degree = 0;

  for (int i = 0; i < NDA; i++) {    // loop over number of distinct abundances
    for (int j = 0; j < f[i]; j++) {  // loop over abundances per class
      for (int nn = 0; nn <= degree; nn++) {
        for (int mm = 1; mm <= g[i]; mm++) {
          if (K[nn] > 0) {
            poly2[nn + mm] += T[i][mm] * K[nn];
          }
        }
        degree += g[i];
        for (int nn = 0; nn <= degree; nn++) {
          K[nn] = (poly2[nn] / powl(10, (4500.0/SPP)));
          poly2[nn] = 0.0;
        }
      }
    }
  }

    // now K[A]=ln(K(D,A)/10^4500) in Etienne's paper
    for (int i = static_cast<int>(SPP); i <= J; i++) {
      K[i] = logl(K[i]);
    }
    // for(i=0;i<NDA;i++) delete[] T[i];
    // delete[] T;
    // delete[] poly2;
    // delete[] f;
    // delete[] g;

  // search of "infinite" values in K[A]
  int borneinf = SPP - 1;
  int bornesup = J + 1;
  long double maxlog = 11333.2;
  bool infinity = false;
  for (int i = SPP; i <= J; i++) {
    if ((K[i] > maxlog) || (K[i] < -maxlog)) {
      infinity = true;
      break;
    }
    borneinf++;
  }    // after that, borneinf=indice next to infinity but before

  for (int i = 0; i <= J - SPP; i++) {
    if ((K[J - i] > maxlog) || (K[static_cast<int>(J) - i] < -maxlog)) {
      infinity = true;
      break;
    }
    bornesup--;
  }    // after that, bornesup=indice next to infinity but after
  if (infinity == true) {
    // cerr << "WARNING : the sample is too large to compute an exact
    // likelihood, the program is thus doing approximations. The following
    // results are to be taken with caution"<< endl;
    // cerr << "Value of A above which K(D,A) is
    // computed approximately ="<< borneinf << endl;
    // cerr << "Value of A below which K(D,A) is
    // computed approximately ="<< bornesup << endl;

    // fitting of the infinite values of K[A] by a polynom of degree 3
    // computing of the derivatives at the critic points

    // Rcpp::Rcout << "Infinity == 1 !! You made it!\n" ;

    if (borneinf > static_cast<int>(K.size())) return K;
    if (borneinf < 1) return K;
    if (bornesup > static_cast<int>((K.size()-1))) return K;
    if (bornesup < 0) return K;
    long double Kprimeinf = K[borneinf] - K[borneinf - 1];


    long double Kprimesup1 = -50.0;   // out of range check
    if ((bornesup + 1) < static_cast<int>(K.size()))
        Kprimesup1 = K[bornesup + 1];
    long double Kprimesup2 = 0.0;    // out of range check
    if (bornesup >= 0.0) Kprimesup2 = K[bornesup];

    // K[bornesup+1]-K[bornesup];
    long double Kprimesup = Kprimesup1 - Kprimesup2;
    // definition of the parameters of the fitted polynom aX^3+bX^2+cX+d
    long double a, b, c, d;
    // inversion of the linear system of equations (with the Gauss method)
    long double borneinf2 = static_cast<long double>(borneinf) *
                            static_cast<long double>(borneinf);
    long double borneinf3 = static_cast<long double>(borneinf2) *
                            static_cast<long double>(borneinf);
    long double bornesup2 = static_cast<long double>(bornesup) *
                            static_cast<long double>(bornesup);
    long double bornesup3 = static_cast<long double>(bornesup2) *
                            static_cast<long double>(bornesup);

    long double ld_borneinf = static_cast<long double>(borneinf);

    d = (Kprimesup - 3 * bornesup2 * K[borneinf]/borneinf3 +
      (2 * bornesup / ld_borneinf - 3 * bornesup2 / borneinf2) *
      (Kprimeinf - 3 * K[borneinf] / ld_borneinf) -
      ((1 + 3 * bornesup2 / borneinf2 - 4 * bornesup/ld_borneinf) /
        (bornesup - 2 * bornesup2 / ld_borneinf +
          bornesup3 / borneinf2)) * (K[bornesup] - bornesup3 *
          K[borneinf] / borneinf3 + (bornesup2 / ld_borneinf -
          bornesup3/borneinf2) * (Kprimeinf - 3 * K[borneinf] /
            ld_borneinf))) / ((6 * bornesup2/borneinf3) -
              (6 * bornesup / borneinf2) - ((1 + 3 * bornesup2/borneinf2 -
              4 * bornesup / ld_borneinf) / (bornesup - 2 * bornesup2 /
                ld_borneinf + bornesup3 / borneinf2)) * (1 - 3 * bornesup2 /
                  borneinf2 + 2 * bornesup3/borneinf3));

    c = ((K[bornesup] - bornesup3 * K[borneinf] / borneinf3 +
      (bornesup2 / ld_borneinf - bornesup3 / borneinf2) *
      (Kprimeinf - 3 * K[borneinf] / ld_borneinf)) - d *
      (1 - 3 * bornesup2/borneinf2 + 2 * bornesup3/borneinf3)) /
        (bornesup - 2 * bornesup2 / ld_borneinf +
          bornesup3 / borneinf2);
    b = (Kprimeinf - 3 * K[borneinf] / ld_borneinf + 2 * c +
      3 * d / ld_borneinf) / (0.0 - ld_borneinf);
    a = (K[borneinf] - b * borneinf2 -
      c * ld_borneinf - d) / borneinf3;

    // reconstruction of K[A] with the fitted polynom
    for (int i = borneinf + 1; i < bornesup; i++) {
      K[i]=(a * i * i * i + b * i * i + c * i + d);
    }
  }

  return K;
}




// [[Rcpp::export]]
NumericVector calcKDA(NumericVector A) {
#ifdef __arm64__
  // long doubles are not supported on ARM / M1 CPU
  // so we have to work around this.
  size_t numspecies = A.size();
  std::vector< size_t > Abund(numspecies);

  for (size_t s = 0; s < numspecies; ++s) {
    if (s > Abund.size()) break;
    Abund[s] = A[s];
  }

  std::vector<double> K = calcLogKDA_arm(numspecies, Abund);
  int sizeofK = K.size();
  Rcpp::NumericVector out(sizeofK);
  for (int i = 0; i < sizeofK; ++i) {
    out[i] = K[i] + 4500.0 * logl(10);
  }

  return out;


#else
  // convert abundances from A to Species
  int numspecies = A.size();
  std::vector<int> Abund(numspecies);  // Abund = new int[numspecies];

  int J = 0;
  for (size_t s = 0; s < numspecies; ++s) {
    if (s > Abund.size()) break;
    Abund[s] = A[s];
    J += Abund[s];
  }
  // call calcLogKDA
  auto K = calcLogKDA(numspecies, Abund);

  // int sizeofK =  J + 1; //I hope!
  NumericVector out(K.size());
  for (size_t i = 0; i < K.size(); ++i) {
    out[i] = K[i] + 4500.0 * logl(10);
  }

  return out;
#endif
}
