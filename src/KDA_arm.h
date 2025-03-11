// Copyright 2015 - 2025 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#ifndef SRC_KDA_ARM_H_
#define SRC_KDA_ARM_H_

#include <cmath>
#include <vector>
#include <algorithm>

#include <Rcpp.h>

class log_val {
 public:
  log_val() {
    is_zero_ = true;
    log_val_ = -1e10;
  }

  explicit log_val(double v) : log_val_(v) {
    is_zero_ = false;
  }

  void reset() {
    is_zero_ = true;
    log_val_ = -1e10;
  }

  void set_with_norm(double v) {
    log_val_ = log(v);
    is_zero_ = false;
  }

  void set_with_log(double l_v) {
    log_val_ = l_v;
    is_zero_ = false;
  }

  bool is_zero() const {
    return is_zero_;
  }

  double get_log_val() const {
    return log_val_;
  }

  log_val& operator+=(const log_val& other) {
    if (is_zero_) {
      if (other.is_zero()) {
        return *this;
      }
      log_val_ = other.get_log_val();
      is_zero_ = false;
      return *this;
    }

    auto l_a = this->get_log_val();
    auto l_b = other.get_log_val();
    log_val_ = l_a + log(1.0 + exp(l_b - l_a));
    return *this;
  }

  log_val& operator+=(const double& other) {
    if (is_zero_) {
      log_val_ = other;
      is_zero_ = false;
      return *this;
    }

    auto l_a = this->get_log_val();
    auto l_b = other;
    log_val_ = l_a + log(1.0 + exp(l_b - l_a));
    return *this;
  }

 private:
  double log_val_;
  bool is_zero_;
};


//// Base of this code by
//   J. Chave and F. Jabot
///  last update 05-23-2008

std::vector<double> calcLogKDA_arm(size_t numspecies,
                                   std::vector<size_t> Abund) {
  if (Abund.size() < 1) return std::vector<double>();

  std::sort(Abund.begin(), Abund.end());

  double J = 0;
  double SPP = numspecies;
  for (size_t s = 0; s < SPP; s++) {
    J += Abund[s];
  }

  double MaxA = Abund.back();   // MaxA = Abund[(int)SPP-1];

  // abundance distribution
  // int *Phi = new int[MaxA+1]; for(s=0;s<=MaxA;s++) Phi[s]=0;
  std::vector< size_t > Phi(MaxA + 1, 0);

  for (size_t s = 0; s < SPP; s++) Phi[ Abund[s] ]++;

  // Number of distinct abundances
  size_t NDA = 0;
  for (size_t s = 0; s <= MaxA; s++) {
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

  std::vector<size_t> f(NDA, 0);   // int *f = new int[NDA];
  std::vector<size_t> g(NDA, 0);   // int *g = new int[NDA];
  size_t i = 0;
  // for(s=0;s<NDA;s++) {f[s]=0;g[s]=0;}
  for (int n = 0; n <= MaxA; n++) {
    if (Phi[n] > 0) {
      f[i] = Phi[n];
      g[i] = n;
      i++;
    }
  }

  std::vector< double > fill;
  std::vector< std::vector< double > > T(NDA, fill);

  std::vector< double > T0(g[0] + 1);
  T[0]    = T0;
  T[0][0] = 0;
  T[0][1] = 1;

  if (g[0] != 1) {
    //  double *lS2 = new double[g[0]+1];
    std::vector<double> lS2(g[0] + 1);
    lS2[0] = 0; lS2[1] = 1;
    for (size_t n = 2; n <= g[0]; n++) {
      std::vector<double> lS1(n + 1);   // double *lS1 = new double[n+1];
      for (size_t im = 0; im <= n-1; im++) {
        lS1[im] = lS2[im];
      }
      lS1[n] = 0;
      for (size_t im = 2; im <= n; im++) {
        lS2[im] = lS1[im] + lS1[im - 1] * (im - 1) / (n - 1);
      }
    }
    for (size_t im = 2; im <= g[0]; im++) {
      T[0][im] = lS2[im];
    }
  }


  for (size_t in = 1; in < i; in++) {
    std::vector<double> Tin(g[in] + 1);   // T[in]= new double[g[in]+1];
    T[in] = Tin;
    T[in][0] = 0;
    T[in][1] = 1;

    std::vector<double> lS2(g[in] + 1);   // double *lS2 = new double[g[in]+1];
    for (size_t im = 0; im <= g[in - 1]; im++) {
      lS2[im] = T[in - 1][im];
    }
    for (size_t n = g[in - 1] + 1; n <= g[in]; n++) {
      std::vector<double> lS1(n + 1);   // double *lS1 = new double[n+1];
      for (size_t im = 0; im <= n - 1; im++) {
        lS1[im] = lS2[im];
      }
      lS1[n] = 0;
      for (size_t im = 2; im <= n; im++) {
        lS2[im] = lS1[im] + lS1[im - 1] * (im - 1) / (n - 1);
      }
    }
    for (size_t im = 2; im <= g[in]; im++) {
      T[in][im] = lS2[im];
    }
  }
  // After this stage we have stored in T[i][m] T(g[i],m)
  // with T(n,m) = S(n,m)*S(m,1)/S(n,1) for i>0

  // precompute log values
  for (auto& i : T) {
    for (auto& j : i) {
      j = log(j);
    }
  }

  // cerr << "Start computing ln(K(D,A)) ...\n";
  // SECOND STAGE: compute the K(D,A)
  // I follow Etienne's route. Compute the product of polynomials
  // of length J

  std::vector< log_val > K_prime(J + 2);
  std::vector< log_val > log_poly2(J + 1);

  K_prime[0].set_with_norm(1);

  size_t degree = 0;
  double minus = 4500.0 * log(10) / SPP;

  // loop over number of distinct abundances
  for (size_t i = 0; i < NDA; i++) {
    for (size_t j = 0; j < f[i]; j++) {   // loop over abundances per class
      for (size_t nn = 0; nn <= degree; nn++) {
        for (size_t mm = 1; mm <= g[i]; mm++) {
          if (!K_prime[nn].is_zero()) {
            log_poly2[nn + mm] += T[i][mm] + K_prime[nn].get_log_val();
          }
        }
      }

      degree += g[i];
      for (size_t nn = 0; nn <= degree; nn++) {
        if (!log_poly2[nn].is_zero()) {
          double add_new = log_poly2[nn].get_log_val() - minus;
          K_prime[nn].set_with_log(add_new);
        } else {
          K_prime[nn].reset();
        }

        log_poly2[nn].reset();
      }
    }
  }


  std::vector<double> K(J + 2, 0.0);

  size_t cnt = 0;
  for (const auto& i : K_prime) {
    K[cnt] = static_cast<double>(i.get_log_val());
    if (K[cnt] == -1e10) K[cnt] = 0.0;
    cnt++;
  }

  double borneinf = SPP - 1;
  double bornesup = J + 1;
  double maxlog = 11333.2;
  bool infinity = false;
  for (int i = SPP; i <= J; i++) {
    if (i >= static_cast<int>(K.size())) {
      Rcpp::stop("i > K.size()");;
    }

    if ((K[i] > maxlog) || (K[i] < -maxlog)) {
      infinity = true;
      break;
    }
    borneinf++;
  }    // after that, borneinf=indice next to infinity but before

  for (int i = 0; i <= J - SPP; i++) {
    if (i >= static_cast<int>(K.size())) {
      Rcpp::stop("i > K.size()");
    }

    if ((K[J - i] > maxlog) || (K[static_cast<int>(J) - i] < -maxlog)) {
      infinity = true;
      break;
    }
    bornesup--;
  }    // after that, bornesup=indice next to infinity but after

  if ( infinity == true ) {
    // cerr << "WARNING : the sample is too large to compute an exact
    // likelihood, the program is thus doing approximations. The following
    // results are to be taken with caution"<<endl;
    // cerr << "Value of A above which K(D,A) is computed approximately ="
    // <<borneinf<<endl;
    // cerr << "Value of A below which K(D,A) is computed approximately ="
    // <<bornesup<<endl;

    // fitting of the infinite values of K[A] by a polynom of degree 3
    // computing of the derivatives at the critic points

    if (static_cast<int>(borneinf) >= static_cast<int>(K.size())) {
      Rcpp::stop("borneinf > K.size()");
    }
    if (bornesup >= static_cast<int>((K.size()))) {
      Rcpp::stop("bornesup > K.size()");
    }

    int index1 = static_cast<int>(borneinf);
    int index2 = static_cast<int>(borneinf - 1);

    if (index1 >= static_cast<int>(K.size()) || index1 < 0) {
      Rcpp::stop("index1 out of range");
    }
    if (index2 >= static_cast<int>(K.size()) || index2 < 0) {
      Rcpp::stop("index2 out of range");
    }

    double Kprimeinf = K[index1] - K[index2];


    double Kprimesup1 = -50.0;  // out of range check

    int bornesup_one = static_cast<int>(bornesup) + 1;

    if (bornesup_one < static_cast<int>(K.size()) && bornesup_one >= 0)
      Kprimesup1 = K[bornesup_one];

    double Kprimesup2 = 0.0;   // out of range check

    if (static_cast<int>(bornesup) >= 0)
      Kprimesup2 = K[static_cast<int>(bornesup)];

    double Kprimesup = Kprimesup1 - Kprimesup2;  // K[bornesup+1]-K[bornesup];
    // definition of the parameters of the fitted polynom aX^3+bX^2+cX+d
    double a, b, c, d;
    // inversion of the linear system of equations (with the Gauss method)
    double borneinf2 = static_cast<double>(borneinf) *
                       static_cast<double>(borneinf);
    double borneinf3 = static_cast<double>(borneinf2) *
                       static_cast<double>(borneinf);
    double bornesup2 = static_cast<double>(bornesup) *
                       static_cast<double>(bornesup);
    double bornesup3 = static_cast<double>(bornesup2) *
                       static_cast<double>(bornesup);

    double d_borneinf = static_cast<double>(borneinf);

    d = (Kprimesup - 3 * bornesup2 * K[borneinf]/borneinf3 +
        (2 * bornesup / d_borneinf - 3 * bornesup2 / borneinf2) *
        (Kprimeinf - 3 * K[borneinf] / d_borneinf) -
        ((1 + 3 * bornesup2 / borneinf2 -
        4 * bornesup / d_borneinf) / (bornesup -
        2 * bornesup2 / d_borneinf + bornesup3 / borneinf2)) *
        (K[bornesup] - bornesup3 * K[borneinf] / borneinf3 +
        (bornesup2 / d_borneinf - bornesup3 / borneinf2) *
        (Kprimeinf - 3 * K[borneinf]/d_borneinf))) /
        ((6 * bornesup2/borneinf3) - (6 * bornesup / borneinf2) -
        ((1 + 3 * bornesup2 / borneinf2 - 4 * bornesup / d_borneinf) /
        (bornesup - 2 * bornesup2 / d_borneinf + bornesup3 / borneinf2)) *
        (1 - 3 * bornesup2 / borneinf2 + 2 * bornesup3 / borneinf3));

    c = ((K[bornesup] - bornesup3 * K[borneinf] / borneinf3 +
        (bornesup2 / d_borneinf - bornesup3 / borneinf2) * (Kprimeinf -
        3 * K[borneinf] / d_borneinf)) - d * (1 -
        3 * bornesup2 / borneinf2 + 2 * bornesup3 / borneinf3)) /
        (bornesup - 2 * bornesup2 / d_borneinf + bornesup3 / borneinf2);
    b = (Kprimeinf -
        3 * K[borneinf] /
          d_borneinf +
          2 * c +
          3 * d / d_borneinf) / (0.0 - d_borneinf);
    a = (K[borneinf] -
        b * borneinf2 -
        c * d_borneinf - d) /
        borneinf3;

    // reconstruction of K[A] with the fitted polynom
    for (int i = borneinf + 1; i < bornesup; i++) {
      K[i] = (a * i * i * i + b * i * i + c * i + d);
    }
  }

  return K;
}

#endif  // SRC_KDA_ARM_H_
