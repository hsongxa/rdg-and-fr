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

#include "reference_segment.h"

int test_reference_segment()
{
  using namespace rdg;

  // linear element
  reference_segment<double> rs1(1);
  auto m_matrix = rs1.mass_matrix();
  auto d_matrix = rs1.derivative_matrix_wrt_r();

  std::cout << "M matrix: " << std::endl << m_matrix << std::endl;
  std::cout << "D matrix: " << std::endl << d_matrix << std::endl;

  // test the SBP property
  auto s_matrix = m_matrix * d_matrix;
  auto s_transpose = d_matrix.transpose() * m_matrix;
  auto sbp = s_matrix + s_transpose;

  std::cout << "Summation by parts: " << std::endl << sbp << std::endl;

  // high order element
  reference_segment<double> rs6(6);
  m_matrix = rs6.mass_matrix();
  d_matrix = rs6.derivative_matrix_wrt_r();

  std::cout << "M matrix: " << std::endl << m_matrix << std::endl;
  std::cout << "D matrix: " << std::endl << d_matrix << std::endl;

  // test the SBP property
  s_matrix = m_matrix * d_matrix;
  s_transpose = d_matrix.transpose() * m_matrix;
  sbp = s_matrix + s_transpose;

  std::cout << "Summation by parts: " << std::endl << sbp << std::endl;

  return 0;
}
