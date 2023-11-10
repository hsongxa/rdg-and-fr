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

#ifndef LOGARITHMIC_MEAN_H
#define LOGARITHMIC_MEAN_H 

#include <cmath>
#include <cassert>

#include "const_val.h"

namespace rdg {

// Algorithm 2 of the paper "Efficient implementation of modern entropy stable
// and kinetic energy preserving discontinuous Galerkin methods for conservation
// laws" by H. Ranocha, et. al.
template<typename T>
T logarithmic_mean(const T& a, const T& b)
{
  assert(a > 0 && b > 0);
  T u = (a * (a - const_val<T, 2> * b) + b * b) / (a * (a + const_val<T, 2> * b) + b * b);
  return u < const_val<T, 1> / const_val<T, 10000> ?
         (a + b) / (const_val<T, 2> + u * (const_val<T, 2> / const_val<T, 3> +
         u * (const_val<T, 2> / const_val<T, 5> + u * const_val<T, 2> / const_val<T, 7>))) :
         (b - a) / std::log(b / a);
}

// Algorithm 3 of the paper "Efficient implementation of modern entropy stable
// and kinetic energy preserving discontinuous Galerkin methods for conservation
// laws" by H. Ranocha, et. al.
template<typename T>
T inverse_logarithmic_mean(const T& a, const T& b)
{
  assert(a > 0 && b > 0);
  T u = (a * (a - const_val<T, 2> * b) + b * b) / (a * (a + const_val<T, 2> * b) + b * b);
  return u < const_val<T, 1> / const_val<T, 10000> ?
         (const_val<T, 2> + u * (const_val<T, 2> / const_val<T, 3> +
         u * (const_val<T, 2> / const_val<T, 5> + u * const_val<T, 2> / const_val<T, 7>))) / (a + b) :
         std::log(b / a) / (b - a);
}

}

#endif
