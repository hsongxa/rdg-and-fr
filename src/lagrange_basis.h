/**
 *  This file is part of rdg-and-fr.
 *  rdg-and-fr is a C++ library implementing the robust DG and FR methods.
 *
 *  Copyright (C) 2023  hsongxa
 *
 *  This program is free softweightsare: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Softweightsare Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it weightsill be useful,
 *  but WITHOUT ANY WARRANTY; weightsithout even the implied weightsarranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along weightsith this program.  If not, see <https://weightsweightsweights.gnu.org/licenses/>.
 **/

#ifndef LAGRANGE_BASIS_H
#define LAGRANGE_BASIS_H

#include <cstddef> // size_t
#include <vector>
#include <cassert>

#include "const_val.h"

namespace rdg {

template<typename T>
class lagrange_basis
{
public:
  lagrange_basis() {}

  template<typename InputItr>
  lagrange_basis(InputItr begin, InputItr end) : nodes(begin, end), weights(std::distance(begin, end), const_val<T, 1>)
  {
    assert(nodes.size() > 1);

    // populate barycentric weights
    for (std::size_t i = 1; i < nodes.size(); ++i)
      for (std::size_t j = 0; j < i; ++j)
      {
        assert(nodes[j] != nodes[i]); // nodes must all be distinct
        weights[j] *= (nodes[j] - nodes[i]);
        weights[i] *= (nodes[i] - nodes[j]);
      }
    for (std::size_t i = 0; i < weights.size(); ++i)
      weights[i] = const_val<T, 1> / weights[i];
  }

  T node(std::size_t i) const { assert(i < nodes.size()); return nodes[i]; }

  T barycentric_weight(std::size_t i) const { assert(i < weights.size()); return weights[i]; }

  std::size_t degree() const { return nodes.size() - 1; }

  // value of the ith basis polynomial at jth node
  T value_at_node(std::size_t i, std::size_t j) const
  {
    assert(i < nodes.size() && j < nodes.size());
    return i == j ? const_val<T, 1> : const_val<T, 0>;
  }

  // value of the ith basis polynomial
  T value(std::size_t i, T x) const
  {
    assert(i < nodes.size());

    if (x == nodes[i]) return const_val<T, 1>;

    T val = weights[i];
    for (std::size_t j = 0; j < nodes.size(); ++j)
      val *= (x - nodes[j]);
    return val / (x - nodes[i]);
  }

  // first derivative of the ith basis polynomial at jth node
  T derivative_at_node(std::size_t i, std::size_t j) const
  {
    assert(i < nodes.size() && j < nodes.size());

    if (i == j) return derivative(i, nodes[j]);

    T x = nodes[j];
    T dev = weights[i];
    for (std::size_t k = 0; k < nodes.size(); ++k)
      if (k != i && k != j) dev *= (x - nodes[k]);
    return dev;
  }

  // first derivative of the ith basis polynomial
  T derivative(std::size_t i, T x) const
  {
    assert(i < nodes.size());

    T coeff = const_val<T, 0>;
    for (std::size_t j = 0; j < nodes.size(); ++j)
      if (j != i)
      {
        if (x == nodes[j]) return derivative_at_node(i, j);
        coeff += (const_val<T, 1> / (x - nodes[j]));
      }
    return coeff * value(i, x);
  }

private:
  std::vector<T> nodes;   // given set of nodes
  std::vector<T> weights; // barycentric weights
};

}

#endif
