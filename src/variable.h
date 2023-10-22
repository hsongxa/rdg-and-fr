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

// use boost::tuple not boost::tuple because boost::tuple
// works with boost::zip_iterator (standard library does
// not have zip_iterator yet)
#include <boost/tuple/tuple.hpp>

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
boost::tuple<T, T, T> initialize_variable_to_zero(type<boost::tuple<T, T, T>>)
{ return boost::make_tuple(const_val<T, 0>, const_val<T, 0>, const_val<T, 0>); }

template<typename T>
boost::tuple<T, T, T, T> initialize_variable_to_zero(type<boost::tuple<T, T, T, T>>)
{ return boost::make_tuple(const_val<T, 0>, const_val<T, 0>, const_val<T, 0>, const_val<T, 0>); }

template<typename T>
boost::tuple<T, T, T, T, T> initialize_variable_to_zero(type<boost::tuple<T, T, T, T, T>>)
{ return boost::make_tuple(const_val<T, 0>, const_val<T, 0>, const_val<T, 0>, const_val<T, 0>, const_val<T, 0>); }

// TODO: add more overloads, even for std::tuples...

template<typename T>
T initialize_variable_to_zero() { return initialize_variable_to_zero(type<T>()); }


// operations on tuples

template<typename T1, typename T2, typename T3>
boost::tuple<T1, T2, T3> operator+(const boost::tuple<T1, T2, T3>& v0, const boost::tuple<T1, T2, T3>& v1)
{
  return boost::make_tuple(boost::get<0>(v0) + boost::get<0>(v1),
                           boost::get<1>(v0) + boost::get<1>(v1),
                           boost::get<2>(v0) + boost::get<2>(v1));
}

template<typename T1, typename T2, typename T3, typename T4, typename T5>
boost::tuple<T1, T2, T3, T4, T5> operator+(const boost::tuple<T1, T2, T3, T4, T5>& v0, const boost::tuple<T1, T2, T3, T4, T5>& v1)
{
  return boost::make_tuple(boost::get<0>(v0) + boost::get<0>(v1),
                           boost::get<1>(v0) + boost::get<1>(v1),
                           boost::get<2>(v0) + boost::get<2>(v1),
                           boost::get<3>(v0) + boost::get<3>(v1),
                           boost::get<4>(v0) + boost::get<4>(v1));
}

template<typename T1, typename T2, typename T3>
boost::tuple<T1, T2, T3>& operator+=(boost::tuple<T1, T2, T3>& v0, const boost::tuple<T1, T2, T3>& v1)
{
  boost::get<0>(v0) += boost::get<0>(v1);
  boost::get<1>(v0) += boost::get<1>(v1);
  boost::get<2>(v0) += boost::get<2>(v1);
  return v0;
}

template<typename T1, typename T2, typename T3, typename T4, typename T5>
boost::tuple<T1, T2, T3, T4, T5>& operator+=(boost::tuple<T1, T2, T3, T4, T5>& v0, const boost::tuple<T1, T2, T3, T4, T5>& v1)
{
  boost::get<0>(v0) += boost::get<0>(v1);
  boost::get<1>(v0) += boost::get<1>(v1);
  boost::get<2>(v0) += boost::get<2>(v1);
  boost::get<3>(v0) += boost::get<3>(v1);
  boost::get<4>(v0) += boost::get<4>(v1);
  return v0;
}

template<typename T1, typename T2, typename T3>
boost::tuple<T1, T2, T3> operator-(const boost::tuple<T1, T2, T3>& v)
{
  return boost::make_tuple(-boost::get<0>(v),
                           -boost::get<1>(v),
                           -boost::get<2>(v));
}

template<typename T1, typename T2, typename T3, typename T4, typename T5>
boost::tuple<T1, T2, T3, T4, T5> operator-(const boost::tuple<T1, T2, T3, T4, T5>& v)
{
  return boost::make_tuple(-boost::get<0>(v),
                           -boost::get<1>(v),
                           -boost::get<2>(v),
                           -boost::get<3>(v),
                           -boost::get<4>(v));
}

template<typename T1, typename T2, typename T3>
boost::tuple<T1, T2, T3> operator-(const boost::tuple<T1, T2, T3>& v0, const boost::tuple<T1, T2, T3>& v1)
{
  return boost::make_tuple(boost::get<0>(v0) - boost::get<0>(v1),
                           boost::get<1>(v0) - boost::get<1>(v1),
                           boost::get<2>(v0) - boost::get<2>(v1));
}

template<typename T1, typename T2, typename T3, typename T4, typename T5>
boost::tuple<T1, T2, T3, T4, T5> operator-(const boost::tuple<T1, T2, T3, T4, T5>& v0, const boost::tuple<T1, T2, T3, T4, T5>& v1)
{
  return boost::make_tuple(boost::get<0>(v0) - boost::get<0>(v1),
                           boost::get<1>(v0) - boost::get<1>(v1),
                           boost::get<2>(v0) - boost::get<2>(v1),
                           boost::get<3>(v0) - boost::get<3>(v1),
                           boost::get<4>(v0) - boost::get<4>(v1));
}

template<typename T1, typename T2, typename T3>
boost::tuple<T1, T2, T3>& operator-=(boost::tuple<T1, T2, T3>& v0, const boost::tuple<T1, T2, T3>& v1)
{
  boost::get<0>(v0) -= boost::get<0>(v1);
  boost::get<1>(v0) -= boost::get<1>(v1);
  boost::get<2>(v0) -= boost::get<2>(v1);
  return v0;
}

template<typename T1, typename T2, typename T3, typename T4, typename T5>
boost::tuple<T1, T2, T3, T4, T5>& operator-=(boost::tuple<T1, T2, T3, T4, T5>& v0, const boost::tuple<T1, T2, T3, T4, T5>& v1)
{
  boost::get<0>(v0) -= boost::get<0>(v1);
  boost::get<1>(v0) -= boost::get<1>(v1);
  boost::get<2>(v0) -= boost::get<2>(v1);
  boost::get<3>(v0) -= boost::get<3>(v1);
  boost::get<4>(v0) -= boost::get<4>(v1);
  return v0;
}

template<typename T, typename T1, typename T2, typename T3>
boost::tuple<T1, T2, T3> operator*(T scalar, const boost::tuple<T1, T2, T3>& v)
{
  return boost::make_tuple(scalar * boost::get<0>(v),
                           scalar * boost::get<1>(v),
                           scalar * boost::get<2>(v));
}

template<typename T, typename T1,  typename T2, typename T3, typename T4, typename T5>
boost::tuple<T1, T2, T3, T4, T5> operator*(T scalar, const boost::tuple<T1, T2, T3, T4, T5>& v)
{
  return boost::make_tuple(scalar * boost::get<0>(v),
                           scalar * boost::get<1>(v),
                           scalar * boost::get<2>(v),
                           scalar * boost::get<3>(v),
                           scalar * boost::get<4>(v));
}

template<typename T, typename T1, typename T2, typename T3>
boost::tuple<T1, T2, T3>& operator*=(boost::tuple<T1, T2, T3>& v, T scalar)
{
  boost::get<0>(v) *= scalar;
  boost::get<1>(v) *= scalar;
  boost::get<2>(v) *= scalar;
  return v;
}

template<typename T, typename T1, typename T2, typename T3, typename T4, typename T5>
boost::tuple<T1, T2, T3, T4, T5>& operator*=(boost::tuple<T1, T2, T3, T4, T5>& v, T scalar)
{
  boost::get<0>(v) *= scalar;
  boost::get<1>(v) *= scalar;
  boost::get<2>(v) *= scalar;
  boost::get<3>(v) *= scalar;
  boost::get<4>(v) *= scalar;
  return v;
}

template<typename T, typename T1, typename T2, typename T3>
boost::tuple<T1, T2, T3> operator/(const boost::tuple<T1, T2, T3>& v, T scalar)
{
  return boost::make_tuple(boost::get<0>(v) / scalar,
                           boost::get<1>(v) / scalar,
                           boost::get<2>(v) / scalar);
}

template<typename T, typename T1, typename T2, typename T3, typename T4, typename T5>
boost::tuple<T1, T2, T3, T4, T5> operator/(const boost::tuple<T1, T2, T3, T4, T5>& v, T scalar)
{
  return boost::make_tuple(boost::get<0>(v) / scalar,
                           boost::get<1>(v) / scalar,
                           boost::get<2>(v) / scalar,
                           boost::get<3>(v) / scalar,
                           boost::get<4>(v) / scalar);
}

// TODO: more types or operations...

}

#endif
