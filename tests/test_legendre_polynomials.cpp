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
#include <iterator>

#include "legendre_polynomials.h"

int test_legendre_polynomials()
{
  using namespace rdg;

  constexpr int OrderMax = 40;
  constexpr int N = 501;

  double x[N];
  double delta = 2. / static_cast<double>(N - 1);
  for (int i = 0; i < N; ++i) 
    x[i] = delta * i - 1;

  double val_batch[N * (OrderMax + 1)];
  double dev_batch[N * (OrderMax + 1)];
  std::vector<double> val_batchv;
  std::vector<double> dev_batchv;
  for (int i = 0; i < N; ++i) 
  {
    legendre_polynomial_values(OrderMax, x[i], val_batch + i * (OrderMax + 1));
    legendre_polynomial_derivatives(OrderMax, x[i], dev_batch + i * (OrderMax + 1));

    legendre_polynomial_values(OrderMax, x[i], std::back_inserter(val_batchv));
    legendre_polynomial_derivatives(OrderMax, x[i], std::back_inserter(dev_batchv));
  }

  // make sure array and vector have the same (batch) results
  for (int i = 0; i < N * (OrderMax + 1); ++i)
  {
      if (val_batch[i] != val_batchv[i])
      {
        std::cout << "val_batch = " << val_batch[i] << ", val_batchv = " << val_batchv[i] << std::endl;
        return 1;
      }
      if (dev_batch[i] != dev_batchv[i])
      {
        std::cout << "dev_batch = " << dev_batch[i] << ", dev_batchv = " << dev_batchv[i] << std::endl;
        return 1;
      }
  }

  double val[N];
  double dev[N];
  std::ofstream file;
  for (int order = 0; order <= OrderMax; ++order)
  {
    for (int i = 0; i < N; ++i) 
    {
      val[i] = legendre_polynomial_value(order, x[i]);
      dev[i] = legendre_polynomial_derivative(order, x[i]);

      // verify that values and derivatives are the same as those in the batch
      int index_to_batch = order + i * (OrderMax + 1);
      if (val[i] != val_batch[index_to_batch])
      {
        std::cout << "val = " << val[i] << ", val_batch = " << val_batch[index_to_batch] << std::endl;
        return 1;
      }
      if (dev[i] != dev_batch[index_to_batch])
      {
        std::cout << "dev = " << dev[i] << ", dev_batch = " << dev_batch[index_to_batch] << std::endl;
        return 1;
      }
    }

    // output values to files
    std::string file_name = "PolynomialValues_order_";
    file_name += std::to_string(order);
    file_name += ".txt";
    file.open(file_name);
    file << "#x                  y" << std::endl;
    for (int i = 0; i < N; ++i) 
      file << x[i] << " " << val[i] << std::endl;
    file.close();

    // output derivatives to files
    file_name = "PolynomialDerivatives_order_";
    file_name += std::to_string(order);
    file_name += ".txt";
    file.open(file_name);
    file << "#x                  y" << std::endl;
    for (int i = 0; i < N; ++i) 
      file << x[i] << " " << dev[i] << std::endl;
    file.close();
  }

  return 0;
}
