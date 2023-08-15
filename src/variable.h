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

#ifndef VARIABLE_H
#define VARIABLE_H 

#include <tuple>

#include "const_val.h"

namespace rdg {

// see https://www.fluentcpp.com/2017/08/15/function-templates-partial-specialization-cpp
// this technique, due to Simon Brand, is used to overload return types, i.e., functions
// with no arguments
template<typename T>
struct type {};

template<typename T>
T initialize_variable_to_zero(type<T>) { return const_val<T, 0>; }

template<typename T>
std::tuple<T, T, T> initialize_variable_to_zero(type<std::tuple<T, T, T>>)
{ return std::make_tuple(const_val<T, 0>, const_val<T, 0>, const_val<T, 0>); }

template<typename T>
std::tuple<T, T, T, T> initialize_variable_to_zero(type<std::tuple<T, T, T, T>>)
{ return std::make_tuple(const_val<T, 0>, const_val<T, 0>, const_val<T, 0>, const_val<T, 0>); }

template<typename T>
std::tuple<T, T, T, T, T> initialize_variable_to_zero(type<std::tuple<T, T, T, T, T>>)
{ return std::make_tuple(const_val<T, 0>, const_val<T, 0>, const_val<T, 0>, const_val<T, 0>, const_val<T, 0>); }

// TODO: add more overloads, even for boost_tuples...

template<typename T>
T initialize_variable_to_zero() { return initialize_variable_to_zero(type<T>()); }

}

#endif
