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

#ifndef FLUX_EULER_1D_H
#define FLUX_EULER_1D_H 

#include <cassert>
#include <cmath>

// boost::tuple works with boost::zip_iterator
#include <boost/tuple/tuple.hpp>

#include "const_val.h"

template<typename T>
class flux_euler_1d
{
public:
  using value_type = T;
  using variable_type = boost::tuple<T, T, T>;

  explicit flux_euler_1d(T gamma) : m_gamma(gamma) {}

  variable_type physical_flux(const variable_type& var) const
  {
    T rho, rhou, E;
    boost::tie(rho, rhou, E) = var;
    assert(rho > 0);

    T u = rhou / rho;
    T p = (m_gamma - rdg::const_val<T, 1>) * (E - rhou * u / rdg::const_val<T, 2>);

    return boost::make_tuple(rhou, rhou * u + p, (E + p) * u);
  }

  // see the KG flux in the paper "Split Form Nodal Discontinuous Galerkin
  // Schemes with Summation-By-Parts Property for the Compressible Euler
  // Equations" by G.J. Gassner, A.R. Winters, and D. Kopriva, 2016
  variable_type numerical_volume_flux(const variable_type& var_minus,
                                      const variable_type& var_plus) const
  {
    T rho_minus, rhou_minus, E_minus, rho_plus, rhou_plus, E_plus;
    boost::tie(rho_minus, rhou_minus, E_minus) = var_minus;
    boost::tie(rho_plus, rhou_plus, E_plus) = var_plus;
    assert(rho_minus > 0 && rho_plus > 0);

    // averages
    T rho = (rho_minus + rho_plus) / rdg::const_val<T, 2>;

    T u_minus = rhou_minus / rho_minus;
    T u_plus = rhou_plus / rho_plus;
    T u = (u_minus + u_plus) / rdg::const_val<T, 2>;

    T p_minus = (m_gamma - rdg::const_val<T, 1>) * (E_minus - rhou_minus * u_minus / rdg::const_val<T, 2>);
    T p_plus = (m_gamma - rdg::const_val<T, 1>) * (E_plus - rhou_plus * u_plus / rdg::const_val<T, 2>);
    T p = (p_minus + p_plus) / rdg::const_val<T, 2>;

    T e = (E_minus / rho_minus + E_plus / rho_plus) / rdg::const_val<T, 2>;

    return boost::make_tuple(rho * u, rho * u * u + p, (rho * e + p) * u);
  }

  // symmetric part plus stabilization part
  // NOTE: unit normal vectors of the faces of a 1d element
  // NOTE: degenerate to a sign, i.e., -1 or +1
  variable_type numerical_surface_flux(const variable_type& var_minus,
                                       const variable_type& var_plus, T sign_minus) const
  {
    T rho_minus, rhou_minus, E_minus, rho_plus, rhou_plus, E_plus;
    boost::tie(rho_minus, rhou_minus, E_minus) = var_minus;
    boost::tie(rho_plus, rhou_plus, E_plus) = var_plus;
    assert(rho_minus > 0 && rho_plus > 0);

    T u_minus = rhou_minus / rho_minus;
    T p_minus = (m_gamma - rdg::const_val<T, 1>) * (E_minus - rhou_minus * u_minus / rdg::const_val<T, 2>); 
    T LF_minus = std::abs(u_minus) + std::sqrt(m_gamma * p_minus / rho_minus);

    T u_plus = rhou_plus / rho_plus;
    T p_plus = (m_gamma - rdg::const_val<T, 1>) * (E_plus - rhou_plus * u_plus / rdg::const_val<T, 2>); 
    T LF_plus = std::abs(u_plus) + std::sqrt(m_gamma * p_plus / rho_plus);

    // averages - symmetric part
    T rho = (rho_minus + rho_plus) / rdg::const_val<T, 2>;
    T u = (u_minus + u_plus) / rdg::const_val<T, 2>;
    T p = (p_minus + p_plus) / rdg::const_val<T, 2>;
    T e = (E_minus / rho_minus + E_plus / rho_plus) / rdg::const_val<T, 2>;

    // jumps for the stabilization part
    T jump_rho = sign_minus < 0 ? rho_plus - rho_minus : rho_minus - rho_plus;
    T jump_rhou = sign_minus < 0 ? rhou_plus - rhou_minus : rhou_minus - rhou_plus;
    T jump_E = sign_minus < 0 ? E_plus - E_minus : E_minus - E_plus;

    // local Lax-Friedrichs flux
    T LF = LF_minus > LF_plus ? LF_minus / rdg::const_val<T, 2> : LF_plus / rdg::const_val<T, 2>;
    return boost::make_tuple(rho * u + LF * jump_rho,
                             rho * u * u + p + LF * jump_rhou,
                             (rho * e + p) * u + LF * jump_E);
  }

private:
  T m_gamma;
};

#endif
