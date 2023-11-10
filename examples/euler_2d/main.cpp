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

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp> // boost::tuple works with boost::zip_iterator

#include "euler_2d.h"
#include "explicit_runge_kutta.h"

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

  int numCells = 1024;
  int order = 2;
  if (argc > 1)
  {
    numCells = std::atoi(argv[1]);
    order = std::atoi(argv[2]);
  }

  euler_2d<double> op(numCells, order);

  // node positions and initial conditions
  int numNodes = op.num_nodes();
  std::vector<double> x(numNodes);
  std::vector<double> d(numNodes); // density rho
  std::vector<double> m(numNodes); // momentum rhou
  std::vector<double> e(numNodes); // energy
  auto varItr = boost::make_zip_iterator(boost::make_tuple(d.begin(), m.begin(), e.begin()));
  op.initialize_dofs(x.begin(), varItr);

  // allocate work space for the Runge-Kutta loop
  std::vector<double> d1(numNodes);
  std::vector<double> m1(numNodes);
  std::vector<double> e1(numNodes);
  auto var1Itr = boost::make_zip_iterator(boost::make_tuple(d1.begin(), m1.begin(), e1.begin()));
  std::vector<double> d2(numNodes);
  std::vector<double> m2(numNodes);
  std::vector<double> e2(numNodes);
  auto var2Itr = boost::make_zip_iterator(boost::make_tuple(d2.begin(), m2.begin(), e2.begin()));
  std::vector<double> d3(numNodes);
  std::vector<double> m3(numNodes);
  std::vector<double> e3(numNodes);
  auto var3Itr = boost::make_zip_iterator(boost::make_tuple(d3.begin(), m3.begin(), e3.begin()));
  std::vector<double> d4(numNodes);
  std::vector<double> m4(numNodes);
  std::vector<double> e4(numNodes);
  auto var4Itr = boost::make_zip_iterator(boost::make_tuple(d4.begin(), m4.begin(), e4.begin()));
  std::vector<double> d5(numNodes);
  std::vector<double> m5(numNodes);
  std::vector<double> e5(numNodes);
  auto var5Itr = boost::make_zip_iterator(boost::make_tuple(d5.begin(), m5.begin(), e5.begin()));

  // time advancing loop
  int maxNumTS = 10000;
  double T = 0.2;
  double t = 0.0;
  double dt = op.timestep_size(varItr);
  std::cout << "dt = " << dt << std::endl;

  auto t0 = std::chrono::system_clock::now();
  int numTS = 0;
  while (t < T && numTS < maxNumTS)
  {
    rdg::rk4(varItr, numNodes, t, dt, op, var1Itr, var2Itr, var3Itr, var4Itr, var5Itr);
    t += dt;
    numTS++;

    dt = op.timestep_size(varItr);
    if ((t + dt) > T) dt = T - t;
    std::cout << "t = " << t << ", next dt = " << dt << std::endl;
  }
  auto t1 = std::chrono::system_clock::now();

  // output to visualize
  std::ofstream file;
  file.open("IsentropicVortexProblem.txt");
  file.precision(std::numeric_limits<double>::digits10);
  file << "#         x         rho" << std::endl;
  for(int i = 0; i < numNodes; ++i)
    file << x[i] << "  " << d[i] << std::endl;
  file << std::endl;
  file << "#         x         u" << std::endl;
  for(int i = 0; i < numNodes; ++i)
    file << x[i] << "  " << m[i] / d[i] << std::endl;
  file << std::endl;
  file << "#         x         p" << std::endl;
  for(int i = 0; i < numNodes; ++i)
  {
    double p = (op.gamma() - 1.) * (e[i] - m[i] * m[i] / (2. * d[i]));
    file << x[i] << "  " << p << std::endl;
  }
  file.close();

  return 0;
}
