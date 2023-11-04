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

#ifndef CONVECTIVE_FLUX_DIV_1D
#define CONVECTIVE_FLUX_DIV_1D 

#include <cassert>
#include <vector>

#include "const_val.h"
#include "variable.h"
#include "dense_matrix.h"
#include "reference_segment.h"

namespace rdg {

// element-wise calculations of divergence of convective flux in one dimensional space
//
// NOTE: Different from div_1d which calculates divergence of the input
// NOTE: variable at their locations, flux_div_1d calculates divergence
// NOTE: of the flux function of the input variable at collocations of
// NOTE: the variable, i.e., no over-integration is used thanks to the
// NOTE: robust DG schemes with flux differencing.
template<typename REFE, typename FLUX> // REFE - 1D reference element 
class convective_flux_div_1d                      // FLUX - flux calculators and associated types
{
public:
  using T = typename FLUX::value_type;

  convective_flux_div_1d(const REFE& elem, const FLUX& flux) : m_ref_elem(&elem), m_flux_op(&flux) {}

  template<typename ZipItr, typename FItr, typename Itr>
  void apply(ZipItr ins, FItr surf_fluxes, T J, Itr outs);

private:
  using V = typename FLUX::variable_type;

  const REFE*     m_ref_elem;
  const FLUX*     m_flux_op;
};

// NOTE: the implementation for 1D is different from 2D & 3D in the following:
// NOTE: 1) the contravariant basis is a constant scalar and cancels with J (their product is one);
// NOTE: 2) the face nodes are hard coded to be 0 and num_nodes() - 1; and
// NOTE: 3) the face mass matrix degenerates to scalar 1
template<typename REFE, typename FLUX> template<typename ZipItr, typename FItr, typename Itr>
void convective_flux_div_1d<REFE, FLUX>::apply(ZipItr ins, FItr surf_fluxes, T J, Itr outs)
{
  assert(J > 0);

  std::size_t N = m_ref_elem->num_nodes();
  std::vector<V> vol_fluxes(N * N);

  // volume integration
  // NOTE: numerical volume fluxes must be consistent and symmetric
  auto D = m_ref_elem->derivative_matrix_wrt_r();
  for(std::size_t i = 0; i < N; ++i)
  {
    for(std::size_t j = 0; j < i; ++j)
      vol_fluxes[i * N + j] = vol_fluxes[j * N + i];

    vol_fluxes[i * N + i] = m_flux_op->physical_flux(*(ins + i));

    for(std::size_t j = i + 1; j < N; ++j)
      vol_fluxes[i * N + j] = m_flux_op->numerical_volume_flux(*(ins + i), *(ins + j));

    *(outs + i) = initialize_variable_to_zero<V>();
    for(std::size_t j = 0; j < N; ++j)
      *(outs + i) += const_val<T, 2> * D(i, j) * vol_fluxes[i * N + j];
  }

  // plus surface integration lifting
  auto M = m_ref_elem->mass_matrix();
  *outs -= (*surf_fluxes - vol_fluxes[0]) / M(0, 0);
  surf_fluxes++;
  *(outs + N - 1) -= (vol_fluxes[N * N - 1] - *surf_fluxes) / M(N - 1, N - 1);

  // divide by J
  T invJ = const_val<T, 1> / J;
  for(std::size_t i = 0; i < N; ++i) *(outs + i) *= invJ;
}

}

#endif
