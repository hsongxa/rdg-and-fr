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
#include <fstream>
#include <string>
#include <vector>

#include "lagrange_basis.h"

int test_lagrange_basis()
{
  using namespace rdg;

  double nodes[]{-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
  lagrange_basis<double> basis(nodes, nodes + 7);

  // check values and derivatives
  std::vector<double> vals;
  std::vector<double> vals_at_node;
  std::vector<double> devs;
  std::vector<double> devs_at_node;

  int degree = basis.degree();
  for (int i = 0; i <= degree; ++i)
    for (int j = 0; j <= degree; ++j)
    {
      vals.push_back(basis.value(i, nodes[j]));
      vals_at_node.push_back(basis.value_at_node(i, j));
      if (vals[vals.size() - 1] != vals_at_node[vals_at_node.size() - 1])
      {
        std::cout << "basis = " << i << ", node = " << j << " found inconsistent values!" << std::endl;
        return 1;
      }

      devs.push_back(basis.derivative(i, nodes[j]));
      devs_at_node.push_back(basis.derivative_at_node(i, j));
      if (devs[devs.size() - 1] != devs_at_node[devs_at_node.size() - 1])
      {
        std::cout << "basis = " << i << ", node = " << j << " found inconsistent derivatives!" << std::endl;
        return 1;
      }
    }

  constexpr int N = 501;
  double x[N];
  for (int p = 0; p < N; ++p) 
    x[p] = p * 2. / static_cast<double>(N - 1) - 1.;

  double val[N];
  double dev[N];
  std::ofstream file;
  for (int b = 0; b <= degree; ++b)
  {
    for (int p = 0; p < N; ++p) 
    {
      val[p] = basis.value(b, x[p]);
      dev[p] = basis.derivative(b, x[p]);
    }

    // output values to files
    std::string file_name = "PolynomialValues_basis_";
    file_name += std::to_string(b);
    file_name += ".txt";
    file.open(file_name);
    file << "#x                  y" << std::endl;
    for (int p = 0; p < N; ++p) 
      file << x[p] << " " << val[p] << std::endl;
    file.close();

    // output derivatives to files
    file_name = "PolynomialDerivatives_basis_";
    file_name += std::to_string(b);
    file_name += ".txt";
    file.open(file_name);
    file << "#x                  y" << std::endl;
    for (int p = 0; p < N; ++p) 
      file << x[p] << " " << dev[p] << std::endl;
    file.close();
  }

  return 0;
}
