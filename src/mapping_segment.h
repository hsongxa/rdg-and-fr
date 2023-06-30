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

#ifndef MAPPING_SEGMENT
#define MAPPING_SEGMENT

#include <cassert>

#include "const_val.h"

namespace rdg {

template<typename T>
class mapping_segment
{
public:
  static T x_to_r(T x0, T x1, T x)
  { assert(x0 < x1); return (const_val<T, 2> * x - x0 - x1) / (x1 - x0); }

  static T r_to_x(T x0, T x1, T r)
  { return ((const_val<T, 1> - r) * x0 + (const_val<T, 1> + r) * x1) / const_val<T, 2>; }

  static T J(T x0, T x1) { return (x1 - x0) / const_val<T, 2>; }

  static T contravariant_basis(T x0, T x1)
  { assert(x0 < x1); return const_val<T, 2> / (x1 - x0); }
};

}

#endif

