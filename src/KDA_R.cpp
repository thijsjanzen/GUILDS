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

// [[Rcpp::export]]
NumericVector calcKDA(NumericVector A) {
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
}
