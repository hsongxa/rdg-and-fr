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

#ifndef REFERENCE_SEGMENT_H
#define REFERENCE_SEGMENT_H 

#include <cstddef> // size_t
#include <iterator>
#include <cassert>

#include "dense_matrix.h"
#include "lagrange_basis.h"
#include "gauss_lobatto_quadrature.h"

namespace rdg {

template<typename T>
class reference_segment
{
public:
  using matrix_type = dense_matrix<T, false>; // row major

  explicit reference_segment(std::size_t order)
  {
    assert(order > 0);

    std::vector<T> nodes;
    gauss_lobatto_quadrature(order + 1, std::back_inserter(nodes), std::back_inserter(weights));
    basis.set_nodes(nodes.begin(), nodes.end());
  }

  std::size_t num_nodes() const { return basis.num_nodes(); }

  T node_position(std::size_t i) const { return basis.node(i); }

  template<typename OutputItr>
  void node_positions(OutputItr it) const;

  matrix_type mass_matrix() const;

  matrix_type derivative_matrix_wrt_r() const;

private:
  lagrange_basis<T> basis;
  std::vector<T> weights; // quadrature weights
};

template<typename T> template<typename OutputItr>
void reference_segment<T>::node_positions(OutputItr it) const
{
  for (std::size_t i = 0; i < num_nodes(); ++i)
    *it++ = node_position(i);
}

template<typename T> 
typename reference_segment<T>::matrix_type reference_segment<T>::mass_matrix() const
{
  std::size_t size = num_nodes();
  matrix_type result = make_zero_matrix<T, false>(size);
  for (std::size_t i = 0; i < size; ++i)
    result(i, i) = weights[i];
  return result;
}

template<typename T> 
typename reference_segment<T>::matrix_type reference_segment<T>::derivative_matrix_wrt_r() const
{
  std::size_t size = num_nodes();
  matrix_type result(size, size);
  for (std::size_t j = 0; j < size; ++j)
    for (std::size_t i = 0; i < size; ++i)
      result(i, j) = basis.derivative_at_node(j, i);
  return result;
}

}

#endif
