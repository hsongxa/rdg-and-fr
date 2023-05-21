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

#ifndef GAUSS_LOBATTO_QUADRATURE_H
#define GAUSS_LOBATTO_QUADRATURE_H

#include <iterator>
#include <memory>

#include <cmath>
#include <cassert>

#include "const_val.h"

namespace rdg {

// see https://stackoverflow.com/questions/29065760/traits-class-to-extract-containers-value-type-from-a-back-insert-iterator
template<class T>
struct output_iterator_traits : std::iterator_traits<T> {};

template< class OutputIt, class T>
struct output_iterator_traits<std::raw_storage_iterator<OutputIt, T>>
: std::iterator<std::output_iterator_tag, T> {};

template<class Container>
struct output_iterator_traits<std::back_insert_iterator<Container>>
: std::iterator<std::output_iterator_tag, typename Container::value_type> {};

template<class Container>
struct output_iterator_traits<std::front_insert_iterator<Container>>
: std::iterator<std::output_iterator_tag, typename Container::value_type> {};

template<class Container>
struct output_iterator_traits<std::insert_iterator<Container>>
: std::iterator<std::output_iterator_tag, typename Container::value_type> {};

template <class T, class charT, class traits>
struct output_iterator_traits<std::ostream_iterator<T, charT, traits>>
: std::iterator<std::output_iterator_tag, T> {};

template <class charT, class traits>
struct output_iterator_traits<std::ostreambuf_iterator<charT, traits>>
: std::iterator<std::output_iterator_tag, charT> {};


// TODO: implement the Golub-Welsch algorithm to cover arbitrary order
template<typename OutputIteratorP, typename OutputIteratorW>
void gauss_lobatto_quadrature(std::size_t npts, OutputIteratorP it_p, OutputIteratorW it_w)
{
  assert(npts >= 2);

  using P = typename output_iterator_traits<OutputIteratorP>::value_type;
  using W = typename output_iterator_traits<OutputIteratorW>::value_type;

  switch (npts)
  {
    case 2:
      *it_p++ = - const_val<P, 1>;
      *it_w++ = const_val<W, 1>;
      *it_p++ = const_val<P, 1>;
      *it_w++ = const_val<W, 1>;
      return;
    case 3:
      *it_p++ = - const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 3>;
      *it_p++ = const_val<P, 0>;
      *it_w++ = const_val<W, 4> / const_val<W, 3>;
      *it_p++ = const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 3>;
      return;
    case 4:
      *it_p++ = - const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 6>;
      *it_p++ = - std::sqrt(const_val<P, 1> / const_val<P, 5>);
      *it_w++ = const_val<W, 5> / const_val<W, 6>;
      *it_p++ = std::sqrt(const_val<P, 1> / const_val<P, 5>);
      *it_w++ = const_val<W, 5> / const_val<W, 6>;
      *it_p++ = const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 6>;
      return;
    case 5:
      *it_p++ = - const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 10>;
      *it_p++ = - std::sqrt(const_val<P, 3> / const_val<P, 7>);
      *it_w++ = const_val<W, 49> / const_val<W, 90>;
      *it_p++ = const_val<P, 0>;
      *it_w++ = const_val<W, 32> / const_val<W, 45>;
      *it_p++ = std::sqrt(const_val<P, 3> / const_val<P, 7>);
      *it_w++ = const_val<W, 49> / const_val<W, 90>;
      *it_p++ = const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 10>;
      return;
    case 6:
      *it_p++ = - const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 15>;
      *it_p++ = - std::sqrt(const_val<P, 1> / const_val<P, 3> + const_val<P, 2> * std::sqrt(const_val<P, 7>) / const_val<P, 21>);
      *it_w++ = (const_val<W, 14> - std::sqrt(const_val<W, 7>)) / const_val<W, 30>;
      *it_p++ = - std::sqrt(const_val<P, 1> / const_val<P, 3> - const_val<P, 2> * std::sqrt(const_val<P, 7>) / const_val<P, 21>);
      *it_w++ = (const_val<W, 14> + std::sqrt(const_val<W, 7>)) / const_val<W, 30>;
      *it_p++ = std::sqrt(const_val<P, 1> / const_val<P, 3> - const_val<P, 2> * std::sqrt(const_val<P, 7>) / const_val<P, 21>);
      *it_w++ = (const_val<W, 14> + std::sqrt(const_val<W, 7>)) / const_val<W, 30>;
      *it_p++ = std::sqrt(const_val<P, 1> / const_val<P, 3> + const_val<P, 2> * std::sqrt(const_val<P, 7>) / const_val<P, 21>);
      *it_w++ = (const_val<W, 14> - std::sqrt(const_val<W, 7>)) / const_val<W, 30>;
      *it_p++ = const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 15>;
      return;
    case 7:
      *it_p++ = - const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 21>;
      *it_p++ = - std::sqrt(const_val<P, 5> / const_val<P, 11> + const_val<P, 2> * std::sqrt(const_val<P, 5> / const_val<P, 3>) / const_val<P, 11>);
      *it_w++ = (const_val<W, 124> - const_val<W, 7> * std::sqrt(const_val<W, 15>)) / const_val<W, 350>;
      *it_p++ = - std::sqrt(const_val<P, 5> / const_val<P, 11> - const_val<P, 2> * std::sqrt(const_val<P, 5> / const_val<P, 3>) / const_val<P, 11>);
      *it_w++ = (const_val<W, 124> + const_val<W, 7> * std::sqrt(const_val<W, 15>)) / const_val<W, 350>;
      *it_p++ = const_val<P, 0>;
      *it_w++ = const_val<W, 256> / const_val<W, 525>;
      *it_p++ = std::sqrt(const_val<P, 5> / const_val<P, 11> - const_val<P, 2> * std::sqrt(const_val<P, 5> / const_val<P, 3>) / const_val<P, 11>);
      *it_w++ = (const_val<W, 124> + const_val<W, 7> * std::sqrt(const_val<W, 15>)) / const_val<W, 350>;
      *it_p++ = std::sqrt(const_val<P, 5> / const_val<P, 11> + const_val<P, 2> * std::sqrt(const_val<P, 5> / const_val<P, 3>) / const_val<P, 11>);
      *it_w++ = (const_val<W, 124> - const_val<W, 7> * std::sqrt(const_val<W, 15>)) / const_val<W, 350>;
      *it_p++ = const_val<P, 1>;
      *it_w++ = const_val<W, 1> / const_val<W, 21>;
      return;
    default:
      throw "gauss-lobatto quadrature order 8 and above are not implemented!";
  }
}

}

#endif
