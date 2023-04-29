/**
 *  This file is part of rdg-and-fr.
 *  rdg-and-fr is a C++ library implementing the robust DG and FR methods.
 *
 *  Copyright (C) 2023  hsongxa
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **/

#ifndef LEGENDRE_POLYNOMIALS_H
#define LEGENDRE_POLYNOMIALS_H

#include <cmath>
#include <cstddef>

#include "const_val.h"

namespace rdg {

template<typename T>
T gamma(T arg) { return std::tgamma(arg); }

template<typename T>
T sqrt(T arg) { return std::sqrt(arg); }

template<typename T>
T pow(T base, T exp) { return std::pow(base, exp); }

template<typename T>
T legendre_polynomial_value(T alpha, T beta, std::size_t n, T x)
{
  T a1 = alpha + const_val<T, 1>;
  T b1 = beta + const_val<T, 1>;
  T ab2 = alpha + beta + const_val<T, 2>;

  T prev_prev_val = sqrt(gamma(ab2) / gamma(a1) / gamma(b1) / pow(const_val<T, 2>, a1 + beta));
  if (n == 0) return prev_prev_val;

  T prev_val = (const_val<T, 1> / const_val<T, 2>) * prev_prev_val * sqrt((ab2 + const_val<T, 1>) / a1 / b1) * (ab2 * x + alpha - beta);
  if (n == 1) return prev_val;

  T val;
  T a_i = const_val<T, 2> / ab2 * sqrt((a1 + beta) * a1 * b1 / (ab2 - const_val<T, 1>) / (ab2 + const_val<T, 1>));
  for (std::size_t i = 1; i < n; ++i)
  {
    T iab = const_val<T, 2> * i + alpha + beta;
    T b_i = (alpha + beta) * (beta - alpha) / iab / (iab + const_val<T, 2>);
    T a_n = const_val<T, 2> / (iab + const_val<T, 2>) * sqrt((i + const_val<T, 1>) * (i + a1 + beta) * (i + a1) * (i + b1) / (iab + const_val<T, 1>) / (iab + const_val<T, 3>));
    val = ((x - b_i) * prev_val - a_i * prev_prev_val) / a_n;

    a_i = a_n;
    prev_prev_val = prev_val;
    prev_val = val;
  }

  return val;
}

template<typename T, typename OutputIterator>
void legendre_polynomial_values(T alpha, T beta, std::size_t n, T x, OutputIterator it)
{
  T a1 = alpha + const_val<T, 1>;
  T b1 = beta + const_val<T, 1>;
  T ab2 = alpha + beta + const_val<T, 2>;

  T prev_prev_val = sqrt(gamma(ab2) / gamma(a1) / gamma(b1) / pow(const_val<T, 2>, a1 + beta));
  it = prev_prev_val;
  if (n == 0) return;

  T prev_val = (const_val<T, 1> / const_val<T, 2>) * prev_prev_val * sqrt((ab2 + const_val<T, 1>) / a1 / b1) * (ab2 * x + alpha - beta);
  it = prev_val;
  if (n == 1) return;

  T a_i = const_val<T, 2> / ab2 * sqrt((a1 + beta) * a1 * b1 / (ab2 - const_val<T, 1>) / (ab2 + const_val<T, 1>));
  for (std::size_t i = 1; i < n; ++i)
  {
    T iab = const_val<T, 2> * i + alpha + beta;
    T b_i = (alpha + beta) * (beta - alpha) / iab / (iab + const_val<T, 2>);
    T a_n = const_val<T, 2> / (iab + const_val<T, 2>) * sqrt((i + const_val<T, 1>) * (i + a1 + beta) * (i + a1) * (i + b1) / (iab + const_val<T, 1>) / (iab + const_val<T, 3>));
    T val = ((x - b_i) * prev_val - a_i * prev_prev_val) / a_n;
    it = val;

    a_i = a_n;
    prev_prev_val = prev_val;
    prev_val = val;
  }
}

template<typename T>
T legendre_polynomial_derivative(T alpha, T beta, std::size_t n, T x)
{
  if (n == 0) return const_val<T, 0>;

  T val = legendre_polynomial_value(alpha + const_val<T, 1>, beta + const_val<T, 1>, n - 1, x);
  return sqrt(n * (n + alpha + beta + const_val<T, 1>)) * val;
}

template<typename T, typename OutputIterator>
void legendre_polynomial_derivatives(T alpha, T beta, std::size_t n, T x, OutputIterator it)
{
  it = const_val<T, 0>;
  if (n == 0) return;
  
  T input_alpha = alpha;
  T input_beta = beta;

  alpha = alpha + const_val<T, 1>;
  beta = beta + const_val<T, 1>;
  n = n - 1;

  // repeat the code for evaluating values at the new (alpha, beta, n)
  T a1 = alpha + const_val<T, 1>;
  T b1 = beta + const_val<T, 1>;
  T ab2 = alpha + beta + const_val<T, 2>;

  T prev_prev_val = sqrt(gamma(ab2) / gamma(a1) / gamma(b1) / pow(const_val<T, 2>, a1 + beta));
  it = sqrt(input_alpha + input_beta + const_val<T, 2>) * prev_prev_val;
  if (n == 0) return;

  T prev_val = (const_val<T, 1> / const_val<T, 2>) * prev_prev_val * sqrt((ab2 + const_val<T, 1>) / a1 / b1) * (ab2 * x + alpha - beta);
  it = sqrt(const_val<T, 2> * (const_val<T, 2> + input_alpha + input_beta + const_val<T, 1>)) * prev_val;
  if (n == 1) return;

  T a_i = const_val<T, 2> / ab2 * sqrt((a1 + beta) * a1 * b1 / (ab2 - const_val<T, 1>) / (ab2 + const_val<T, 1>));
  for (std::size_t i = 1; i < n; ++i)
  {
    T iab = const_val<T, 2> * i + alpha + beta;
    T b_i = (alpha + beta) * (beta - alpha) / iab / (iab + const_val<T, 2>);
    T a_n = const_val<T, 2> / (iab + const_val<T, 2>) * sqrt((i + const_val<T, 1>) * (i + a1 + beta) * (i + a1) * (i + b1) / (iab + const_val<T, 1>) / (iab + const_val<T, 3>));
    T val = ((x - b_i) * prev_val - a_i * prev_prev_val) / a_n;
    it = sqrt((i + 2) * ((i + 2) + input_alpha + input_beta + const_val<T, 1>)) * val;

    a_i = a_n;
    prev_prev_val = prev_val;
    prev_val = val;
  }
}

}

#endif
