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
 
#ifndef ADVECTION_1D_H
#define ADVECTION_1D_H 

#include <vector>
#include <iterator>
#include <algorithm>
#include <math.h>

#include "uniform_cartesian_mesh_1d.h"
#include "mapping_segment.h"
#include "reference_segment.h"
#include "flux_advection_1d.h"
#include "flux_div_1d.h"

// host code of the problem of linear advection equation in one dimensional space
template<typename T>
class advection_1d
{
public:
  advection_1d(std::size_t numCells, int order)
    : m_numCells(numCells), m_order(order), m_mesh((T)(-M_PI), (T)(M_PI), numCells) {}
  ~advection_1d(){}
  
  T wave_speed() const { return s_waveSpeed; }

  T min_elem_size() const { return (T)(2.L) * (T)(M_PI) / m_numCells; }

  int num_dofs() const { return m_numCells * (m_order + 1); }

  // the layout of DOFs in memory are different for CPU execution and GPU execution;
  // the first iterator sets the DOF positions and the second iterator sets the
  // initial values of the DOFs
  template<typename OutputIterator1, typename OutputIterator2>
  void initialize_dofs(OutputIterator1 it1, OutputIterator2 it2) const;

  // the layout of DOFs in memory are different for CPU execution and GPU execution
  template<typename OutputIterator>
  void exact_solution(T t, OutputIterator it) const;

  // CPU execution of the spatial discrete operator
  template<typename ConstItr, typename Itr>
  void operator()(ConstItr in_cbegin, std::size_t size, T t, Itr out_begin) const;

private:
  template<typename ConstItr>
  void numerical_fluxes(ConstItr cbegin, T t) const; // time t is used for boundary conditions

private:
  using mesh_type         = rdg::uniform_cartesian_mesh_1d<T>;
  using mapping_type      = rdg::mapping_segment;
  using reference_element = rdg::reference_segment<T>;
  using flux_calculator   = rdg::flux_advection_1d<T>;

  // numerical scheme data (could be constants if never change)
  std::size_t m_numCells;
  int m_order;

  // problem definitions
  const mesh_type m_mesh;
  const T s_waveSpeed = (T)(2.L) * (T)(M_PI);

  // work space for numerical fluxes
  mutable std::vector<T> m_numericalFluxes;
};

template<typename T> template<typename OutputIterator1, typename OutputIterator2>
void advection_1d<T>::initialize_dofs(OutputIterator1 it1, OutputIterator2 it2) const
{
  reference_element refElem(m_order);
  std::vector<T> pos;
  refElem.node_positions(std::back_inserter(pos));

  for (std::size_t i = 0; i < m_numCells; ++i)
  {
    auto cell= m_mesh.get_cell(i);
    for (std::size_t j = 0; j < pos.size(); ++j)
    {
      T x = mapping_type::r_to_x(std::get<0>(cell), std::get<1>(cell), pos[j]);
      *it1++ = x;
      *it2++ = std::sin(x);
    }
  }
}

template<typename T> template<typename OutputIterator>
void advection_1d<T>::exact_solution(T t, OutputIterator it) const
{
  reference_element refElem(m_order);
  std::vector<T> pos;
  refElem.node_positions(std::back_inserter(pos));

  for (std::size_t i = 0; i < m_numCells; ++i)
  {
    auto cell= m_mesh.get_cell(i);
    for (std::size_t j = 0; j < pos.size(); ++j)
    {
      T x = mapping_type::r_to_x(std::get<0>(cell), std::get<1>(cell), pos[j]);
      *it++ = std::sin(x - s_waveSpeed * t);
    }
  }
}

template<typename T> template<typename ConstItr>
void advection_1d<T>::numerical_fluxes(ConstItr cbegin, T t) const
{
  flux_calculator fluxCalculator(s_waveSpeed);

  std::size_t numFluxes = m_numCells + 1;
  if (m_numericalFluxes.size() < numFluxes) m_numericalFluxes.resize(numFluxes);

  T a, b;
  int np = m_order + 1; // d.o.f. management !
  for (std::size_t i = 0; i < numFluxes; ++i)
  {
    if (i > 0) a = *(cbegin + (i * np - 1));
    else a = - sin(2.0 * M_PI * t); // inflow boundary condition
    if (i < numFluxes - 1) b = *(cbegin + (i * np));
    else b = *(cbegin + (i * np - 1)); // outflow boundary condition - alternatively, may be set to zero ?
    m_numericalFluxes[i] = fluxCalculator.numerical_surface_flux(a, b, 1);
  }
}

template<typename T> template<typename ConstItr, typename Itr>
void advection_1d<T>::operator()(ConstItr in_cbegin, std::size_t size, T t, Itr out_begin) const
{
  numerical_fluxes(in_cbegin, t);

  rdg::flux_div_1d<reference_element, flux_calculator> divOp{reference_element(m_order), flux_calculator(s_waveSpeed)};

  for (std::size_t cell = 0; cell < m_numCells; ++cell)
  {
    auto cellGeom = m_mesh.get_cell(cell);
    T J = mapping_type::J(std::get<0>(cellGeom), std::get<1>(cellGeom));
    divOp.apply(in_cbegin, m_numericalFluxes.cbegin(), J, out_begin);
  }
}

#endif
