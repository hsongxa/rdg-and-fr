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
#include <iostream>
#include <fstream>
#include <limits>
#include <chrono>

#include "advection_1d.h"
#include "explicit_runge_kutta.h"

double compute_error_norm(double* ref_v, double* v, int size)
{
  double err = 0.0;
  for(int i = 0; i < size; ++i)
    err += (ref_v[i] - v[i]) * (ref_v[i] - v[i]);
  return err / size;
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

  int numCells = 1024;
  int order = 1;
  if (argc > 1)
  {
    numCells = std::atoi(argv[1]);
    order = std::atoi(argv[2]);
  }

  advection_1d<double> op(numCells, order);

  // DOF positions and initial conditions
  int numDOFs = op.num_dofs();
  std::vector<double> x(numDOFs);
  std::vector<double> v(numDOFs);
  op.initialize_dofs(x.begin(), v.begin());

  // allocate work space for the Runge-Kutta loop and the reference solution
  std::vector<double> v1(numDOFs);
  std::vector<double> v2(numDOFs);
  std::vector<double> v3(numDOFs);
  std::vector<double> v4(numDOFs);
  std::vector<double> v5(numDOFs);
  std::vector<double> ref_v(numDOFs);
  
  // time advancing loop
  int totalTSs = 10000;
  double t = 0.0;
  double dt = 0.25 / order / order * op.min_elem_size() / op.wave_speed();

  auto t0 = std::chrono::system_clock::now();
  for (int i = 0; i < totalTSs; ++i)
  {
    rdg::rk4(v.begin(), numDOFs, t, dt, op, v1.begin(), v2.begin(), v3.begin(), v4.begin(), v5.begin());
    t += dt;
  }
  auto t1 = std::chrono::system_clock::now();

  // exact solution
  op.exact_solution(t, ref_v.begin());

  // output the last error
  double errNorm = compute_error_norm(ref_v.data(), v.data(), numDOFs);
  std::cout << "t = " << t << ", error norm = " << errNorm << std::endl;
  std::cout << "time used: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms" << std::endl;

  // output to visualize
  std::ofstream file;
  file.open("Advection1DDataFile.txt");
  file.precision(std::numeric_limits<double>::digits10);
  file << "#         x         y" << std::endl;
  for(int i = 0; i < numDOFs; ++i)
    file << x[i] << "  " << v[i] << std::endl;
  file << std::endl;
  file << "#         x         reference solution" << std::endl;
  for(int i = 0; i < numDOFs; ++i)
    file << x[i] << " " << ref_v[i] << std::endl;

  return 0;
}
