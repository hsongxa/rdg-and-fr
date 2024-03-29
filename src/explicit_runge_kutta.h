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

// NOTE: De-referencing boost::zip_iterator does not directly give boost::tuple;
// NOTE: instead, it gives a boost internal representation. So we need to explicitly
// NOTE: convert to the specific variable_type/boost::tuple, for which we have
// NOTE: operators like +, -, *, /, +=, *=, ..., etc., overloaded (in variable.h).
template <typename T, typename Itr, typename VART = typename Itr::value_type>
void axpy_n(T a, Itr x_cbegin, std::size_t x_size, Itr y_cbegin, Itr out_begin)
{
  assert(x_cbegin != y_cbegin);
  assert(out_begin != x_cbegin && out_begin != y_cbegin);

  for (std::size_t i = 0; i < x_size; ++i) *out_begin++ = a * VART(*x_cbegin++) + VART(*y_cbegin++);
}

// fourth-order explicit Runge-Kutta scheme
template <typename Itr, typename T, typename DiscreteOp>
void rk4(Itr inout, std::size_t size, T t, T dt, const DiscreteOp& op, Itr wk0, Itr wk1, Itr wk2, Itr wk3, Itr wk4)
{
  T half = const_val<T, 1> / const_val<T, 2>;

  op(inout, size, t, wk1);

  axpy_n<T, Itr, typename DiscreteOp::variable_type>(half * dt, wk1, size, inout, wk0);
  op(wk0, size, t + half * dt, wk2);

  axpy_n<T, Itr, typename DiscreteOp::variable_type>(half * dt, wk2, size, inout, wk0);
  op(wk0, size, t + half * dt, wk3);

  axpy_n<T, Itr, typename DiscreteOp::variable_type>(dt, wk3, size, inout, wk0);
  op(wk0, size, t + dt, wk4);

  axpy_n<T, Itr, typename DiscreteOp::variable_type>(const_val<T, 2>, wk2, size, wk1, wk0);
  axpy_n<T, Itr, typename DiscreteOp::variable_type>(const_val<T, 2>, wk3, size, wk4, wk1);
  axpy_n<T, Itr, typename DiscreteOp::variable_type>(dt / const_val<T, 6>, wk0, size, inout, wk2);
  axpy_n<T, Itr, typename DiscreteOp::variable_type>(dt / const_val<T, 6>, wk1, size, wk2, inout);
}

}

#endif
