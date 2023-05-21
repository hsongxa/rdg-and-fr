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

#include <vector>
#include <iterator>
#include <iostream>

#include "gauss_lobatto_quadrature.h"

int test_gauss_lobatto_quadrature()
{
  using namespace rdg;

  std::vector<double> points;
  std::vector<double> weights;

  // Gauss-Lobatto quadrature
  for(int np = 2; np < 8; ++np)
  {
    points.clear();
    weights.clear();
    gauss_lobatto_quadrature(np, std::back_inserter(points), std::back_inserter(weights));

    std::cout << "Gauss-Lobatto quadrature of " << np << " points:" << std::endl;
    for (std::size_t j = 0; j < points.size(); ++j)
      std::cout << "p = " << points[j] << ", w = " << weights[j] << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
