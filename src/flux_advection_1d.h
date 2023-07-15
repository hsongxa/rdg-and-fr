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

#ifndef FLUX_ADVECTION_1D_H
#define FLUX_ADVECTION_1D_H 

#include <cmath>

#include "const_val.h"

namespace rdg {

template<typename T>
class flux_advection_1d
{
public:
  using value_type = T;
  using variable_type = T; // scalar variable

  explicit flux_advection_1d(T velocity) : m_velocity(velocity) {}

  T physical_flux(const T& u) const { return m_velocity * u; }

  // symmetric
  T numerical_volume_flux(const T& u_a, const T& u_b) const
  { return m_velocity * (u_a + u_b) / const_val<T, 2>; }

  // symmetric part plus stabilization part: upwind
  // NOTE: unit normal vectors of the faces of a 1d element
  // NOTE: degenerate to a sign, i.e., -1 or +1
  T numerical_surface_flux(const T& u_a, const T& u_b, T sign_a) const
  {
    T jump = sign_a < 0 ? u_b - u_a : u_a - u_b;
    return numerical_volume_flux(u_a, u_b) + std::fabs(m_velocity) * jump / const_val<T, 2>;
  }

private:
  T m_velocity;
};

}

#endif
