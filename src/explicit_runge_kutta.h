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

#ifndef EXPLICIT_RUNGE_KUTTA_H
#define EXPLICIT_RUNGE_KUTTA_H

#include <cassert>
#include <cstddef>

#include "const_val.h"

namespace rdg {

// axpy
template <typename T, typename ConstItr, typename Itr>
void axpy_n(T a, ConstItr x_cbegin, std::size_t x_size, ConstItr y_cbegin, Itr out_begin)
{
  assert(x_cbegin != y_cbegin);
  assert(out_begin != x_cbegin && out_begin != y_cbegin);

  for (std::size_t i = 0; i < x_size; ++i) *out_begin++ = a * (*x_cbegin++) + (*y_cbegin++);
}

// fourth-order explicit Runge-Kutta scheme
template <typename Itr, typename T, typename DiscreteOp, typename Axpy>
void rk4(Itr inout, std::size_t size, T t, T dt, const DiscreteOp& op, const Axpy& axpy, Itr wk0, Itr wk1, Itr wk2, Itr wk3, Itr wk4)
{
  T half = const_val<T, 1> / const_val<T, 2>;

  op(inout, size, t, wk1);

  axpy(half * dt, wk1, size, inout, wk0);
  op(wk0, size, t + half * dt, wk2);

  axpy(half * dt, wk2, size, inout, wk0);
  op(wk0, size, t + half * dt, wk3);

  axpy(dt, wk3, size, inout, wk0);
  op(wk0, size, t + dt, wk4);

  axpy(const_val<T, 2>, wk2, size, wk1, wk0);
  axpy(const_val<T, 2>, wk3, size, wk4, wk1);
  axpy(dt / const_val<T, 6>, wk0, size, inout, wk2);
  axpy(dt / const_val<T, 6>, wk1, size, wk2, inout);
}

}

#endif
