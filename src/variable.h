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

#include <tuple> // TODO: change to std::tuple when zip_iterators become standard

#include "const_val.h"

namespace rdg {

template<typename... Types>
using tuple = std::tuple<Types...>; // TODO: change to std::tuple when zip_iterators become standard


template<typename T>
T make_zero_variable() { return const_val<T, 0>; }

template<typename tuple<T, T, T>>
tuple<T, T, T> make_zero_variable()
{ return std::make_tuple(const_val<T, 0>, const_val<T, 0>, const_val<T, 0>); }

template<typename tuple<T, T, T, T>>
tuple<T, T, T, T> make_zero_variable()
{ return std::make_tuple(const_val<T, 0>, const_val<T, 0>, const_val<T, 0>, const_val<T, 0>); }

template<typename tuple<T, T, T, T, T>>
tuple<T, T, T, T, T> make_zero_variable()
{ return std::make_tuple(const_val<T, 0>, const_val<T, 0>, const_val<T, 0>, const_val<T, 0>, const_val<T, 0>); }

// TODO: add overloads for other possible variable types here...

}

#endif
