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

#ifndef UNIFORM_CARTESIAN_MESH_H
#define UNIFORM_CARTESIAN_MESH_H 

#include <cassert>
#include <cstddef>
#include <tuple>

namespace rdg {

// a simple implementation of 1D mesh for testing purposes
template<typename T>
class uniform_cartesian_mesh_1d
{
public:
  using point_type = T;

  uniform_cartesian_mesh_1d(T x0, T x1, std::size_t n) : m_x0(x0), m_n(n)
  {
    assert(x0 < x1 && n > 0);
    m_delta = (x1 - x0) / static_cast<T>(n);
  }

  std::size_t num_vertices() const { return m_n + 1; }

  std::size_t num_cells() const { return m_n; }

  point_type get_vertex(std::size_t i) const { return m_x0 + i * m_delta; }

  std::tuple<point_type, point_type> get_cell(std::size_t i) const
  { return std::make_tuple(m_x0 + i * m_delta, m_x0 + (i + 1) * m_delta); }

private:
  T m_x0;
  T m_delta;
  std::size_t m_n;
};

}

#endif
