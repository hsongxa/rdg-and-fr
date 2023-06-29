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

#include <iostream>

#include "mapping_segment.h"

int test_mapping_segment()
{
  using namespace rdg;

  const double A = -2.;
  const double B = -1.;
  double x = -1.25;
  double r = mapping_segment<double>::x_to_r(A, B, x);
  std::cout << "in segment [-2, -1], x = -1.25 is mapped to r = " << r << std::endl;

  x = mapping_segment<double>::r_to_x(A, B, 0.5);
  std::cout << "in segment [-2, -1], r = 0.5 is mapped to x = " << x << std::endl;

  double J = mapping_segment<double>::J(A, B);
  std::cout << "J of the segment [-2, -1] = " << J << std::endl;

  std::cout << "contravariant basis of the segment [-2, -1] = ";
  std::cout << mapping_segment<double>::contravariant_basis(A, B) << std::endl;

  return 0;
}
