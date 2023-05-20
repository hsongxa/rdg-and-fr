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

#include <cstddef> // size_t
#include <cassert>

namespace rdg {

template<typename T>
T legendre_polynomial_value(std::size_t n, T x)
{
  assert(x >= static_cast<T>(-1) && x <= static_cast<T>(1));

  T prev_prev_val = static_cast<T>(1);
  if (n == 0) return prev_prev_val;

  T prev_val = x;
  if (n == 1) return prev_val;

  T val;
  for (std::size_t i = 1; i < n; ++i)
  {
    val = (static_cast<T>(2 * i + 1) * prev_val * x - static_cast<T>(i) * prev_prev_val) /
          static_cast<T>(i + 1);

    prev_prev_val = prev_val;
    prev_val = val;
  }

  return val;
}

template<typename T, typename OutputIterator>
void legendre_polynomial_values(std::size_t n, T x, OutputIterator it)
{
  assert(x >= static_cast<T>(-1) && x <= static_cast<T>(1));

  T prev_prev_val = static_cast<T>(1);
  *it++ = prev_prev_val;
  if (n == 0) return;

  T prev_val = x;
  *it++ = prev_val;
  if (n == 1) return;

  for (std::size_t i = 1; i < n; ++i)
  {
    T val = (static_cast<T>(2 * i + 1) * prev_val * x - static_cast<T>(i) * prev_prev_val) /
            static_cast<T>(i + 1);
    *it++ = val;

    prev_prev_val = prev_val;
    prev_val = val;
  }
}

template<typename T>
T legendre_polynomial_derivative(std::size_t n, T x)
{
  assert(x >= static_cast<T>(-1) && x <= static_cast<T>(1));

  if (x == static_cast<T>(1) || x == static_cast<T>(-1))
    return n % 2 == 0 ? x * static_cast<T>(n * (n + 1)) / static_cast<T>(2) :
                        static_cast<T>(n * (n + 1)) / static_cast<T>(2);

  if (n == 0) return static_cast<T>(0);
  if (n == 1) return static_cast<T>(1);


  T prev_prev_val = static_cast<T>(1);
  T prev_val = x;
  for (std::size_t i = 1; i < n; ++i)
  {
    T val = (static_cast<T>(2 * i + 1) * prev_val * x - static_cast<T>(i) * prev_prev_val) /
            static_cast<T>(i + 1);

    if (i + 1 == n) return (val * x - prev_val) * static_cast<T>(i + 1) / (x * x - static_cast<T>(1));

    prev_prev_val = prev_val;
    prev_val = val;
  }
}

template<typename T, typename OutputIterator>
void legendre_polynomial_derivatives(std::size_t n, T x, OutputIterator it)
{
  assert(x >= static_cast<T>(-1) && x <= static_cast<T>(1));

  if (x == static_cast<T>(1) || x == static_cast<T>(-1))
  {
    for (std::size_t i = 0; i <= n; ++i)
      *it++ = i % 2 == 0 ? x * static_cast<T>(i * (i + 1)) / static_cast<T>(2) :
                           static_cast<T>(i * (i + 1)) / static_cast<T>(2);
    return;
  }

  *it++ = static_cast<T>(0);
  if (n == 0) return;

  *it++ = static_cast<T>(1);
  if (n == 1) return;

  T prev_prev_val = static_cast<T>(1);
  T prev_val = x;
  for (std::size_t i = 1; i < n; ++i)
  {
    T val = (static_cast<T>(2 * i + 1) * prev_val * x - static_cast<T>(i) * prev_prev_val) /
            static_cast<T>(i + 1);

    *it++ = (val * x - prev_val) * static_cast<T>(i + 1) / (x * x - static_cast<T>(1));

    prev_prev_val = prev_val;
    prev_val = val;
  }
}

template<typename T>
T legendre_polynomial_l2_norm(std::size_t n)
{ return static_cast<T>(2) / static_cast<T>(2 * n + 1); }

}

#endif
