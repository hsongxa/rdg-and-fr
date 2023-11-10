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
 
#ifndef EULER_2D_H
#define EULER_2D_H 

#include <vector>
#include <cmath>
#include <cassert>
#include <limits>

// use boost::tuple instead of std::tuple because boost::tuple
// can work with boost::zip_iterator; standard library does not
// have zip_iterator yet
#include <boost/tuple/tuple.hpp>
// operations on tuples
#include "variable.h"

#include "uniform_cartesian_mesh_1d.h"
#include "mapping_segment.h"
#include "reference_segment.h"
#include "flux_euler_2d.h"
#include "convective_flux_div_1d.h"


// host code of the problem of euler equation in one dimensional space
template<typename T>
class euler_2d
{
public:
  euler_2d(std::size_t numCells, int order)
    : m_numCells(numCells), m_order(order), m_mesh((T)(0), (T)(1), numCells) {}
  ~euler_2d(){}

  T gamma() const { return s_gamma; }

  int num_nodes() const { return m_numCells * (m_order + 1); }

  // the layout of DOFs in memory are different for CPU execution and GPU execution;
  // the first iterator sets the node positions and the second iterator sets the
  // initial values of the DOFs
  template<typename OutputIterator1, typename OutputZipIterator2>
  void initialize_dofs(OutputIterator1 it1, OutputZipIterator2 it2) const;

  // suggested next timestep size
  // input is the solution at the current timestep
  template<typename InputZipIterator>
  T timestep_size(InputZipIterator it) const; 

  // CPU execution of the spatial discrete operator
  template<typename ConstZipItr, typename ZipItr>
  void operator()(ConstZipItr in_cbegin, std::size_t size, T t, ZipItr out_begin) const;

  using variable_type = boost::tuple<T, T, T>;

private:
  template<typename ConstItr>
  void numerical_fluxes(ConstItr cbegin, T t) const; // time t is used for boundary conditions

  // re-direct the search for the unary operator- to the rdg name space
  static variable_type negative(const variable_type& v) { return rdg::operator-(v); }

private:
  using mesh_type         = rdg::uniform_cartesian_mesh_1d<T>;
  using mapping_type      = rdg::mapping_segment;
  using reference_element = rdg::reference_segment<T>;
  using flux_calculator   = flux_euler_2d<T>;

  // numerical scheme data (could be constants if never change)
  std::size_t m_numCells;
  int m_order;

  // problem definitions
  const mesh_type m_mesh;
  const T s_gamma = static_cast<T>(1.4);

  // work space for numerical fluxes
  mutable std::vector<variable_type> m_numericalFluxes;
};

template<typename T> template<typename OutputIterator1, typename OutputZipIterator2>
void euler_2d<T>::initialize_dofs(OutputIterator1 it1, OutputZipIterator2 it2) const
{
  reference_element refElem(m_order);
  std::vector<T> pos;
  refElem.node_positions(std::back_inserter(pos));

  // conserved variables, not primary variables
  for (std::size_t i = 0; i < m_numCells; ++i)
  {
    auto cell= m_mesh.get_cell(i);
    for (std::size_t j = 0; j < pos.size(); ++j)
    {
      T x = mapping_type::r_to_x(std::get<0>(cell), std::get<1>(cell), pos[j]);
      *it1++ = x;
      *it2++ = x < 0.5 ?
               boost::make_tuple(static_cast<T>(1), static_cast<T>(0), static_cast<T>(1) / (s_gamma - static_cast<T>(1))) :
               boost::make_tuple(static_cast<T>(0.125), static_cast<T>(0), static_cast<T>(0.1) / (s_gamma - static_cast<T>(1)));
    }
  }
}

template<typename T> template<typename InputZipIterator>
T euler_2d<T>::timestep_size(InputZipIterator it) const
{
  T rho, rhou, E;

  T maxV = std::numeric_limits<T>::lowest();
  for (std::size_t i = 0; i < num_nodes(); ++i)
  {
    boost::tie(rho, rhou, E) = *it++;
    T u = rhou / rho;
    T p = (E - rhou * u / static_cast<T>(2)) * (s_gamma - static_cast<T>(1)); 
    T v = std::abs(u) + std::sqrt(s_gamma * p / rho);
    if (v > maxV) maxV = v;
  }

  return static_cast<T>(0.25) / (maxV * static_cast<T>(m_mesh.num_cells())) / static_cast<T>(m_order);
}

template<typename T> template<typename ConstZipItr>
void euler_2d<T>::numerical_fluxes(ConstZipItr cbegin, T t) const
{
  flux_calculator fluxCalculator(s_gamma);

  std::size_t numFluxes = m_numCells + 1;
  if (m_numericalFluxes.size() < numFluxes) m_numericalFluxes.resize(numFluxes);

  variable_type a, b;
  int np = reference_element(m_order).num_nodes();
  for (std::size_t i = 0; i < numFluxes; ++i)
  {
    if (i > 0) a = *(cbegin + (i * np - 1));
    // inflow boundary condition
    else a = boost::make_tuple(static_cast<T>(1), static_cast<T>(0), static_cast<T>(1) / (s_gamma - static_cast<T>(1)));
    if (i < numFluxes - 1) b = *(cbegin + (i * np));
    // outflow boundary condition
    else b = boost::make_tuple(T(0.125), static_cast<T>(0), T(0.1) / (s_gamma - static_cast<T>(1))); //*(cbegin + (i * np - 1));
    m_numericalFluxes[i] = fluxCalculator.numerical_surface_flux(a, b, 1);
  }
}

template<typename T> template<typename ConstZipItr, typename ZipItr>
void euler_2d<T>::operator()(ConstZipItr in_cbegin, std::size_t size, T t, ZipItr out_begin) const
{
  numerical_fluxes(in_cbegin, t);

  reference_element refElem(m_order);
  flux_calculator fluxCalculator(s_gamma);
  rdg::convective_flux_div_1d<reference_element, flux_calculator> divOp(refElem, fluxCalculator);

  int np = refElem.num_nodes();
  std::vector<variable_type> cellOut(np);
  for (std::size_t cell = 0; cell < m_numCells; ++cell)
  {
    auto cellGeom = m_mesh.get_cell(cell);
    T J = mapping_type::J(std::get<0>(cellGeom), std::get<1>(cellGeom));
    divOp.apply(in_cbegin + np * cell, m_numericalFluxes.cbegin() + cell, J, cellOut.begin());

    for(int i = 0; i < np; ++i) *(out_begin + np * cell + i) = negative(cellOut[i]);
  }
}

#endif
