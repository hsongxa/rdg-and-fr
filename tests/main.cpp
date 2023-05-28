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

#include "unittests.h"

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

//  if (test_legendre_polynomials())
//    std::cout << "test_legendre_polynomials FAILED!!!" << std::endl;
//
//  if (test_gauss_lobatto_quadrature())
//    std::cout << "test_gauss_lobatto_quadrature FAILED!!!" << std::endl;

    if (test_lagrange_basis())
      std::cout << "test_lagrange_basis FAILED!!!" << std::endl;

  // finish
  return 0;
}
