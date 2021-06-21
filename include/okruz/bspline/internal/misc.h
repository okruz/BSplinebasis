#ifndef OKRUZ_BSPLINE_MISC_H
#define OKRUZ_BSPLINE_MISC_H

/*
 * ########################################################################
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
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ########################################################################
 */

/*!
 * The methods in the internal namespace are internal helper method for this
 * library.
 */

namespace okruz::bspline::internal {

/*!
 * Creates an std::array<T, size> with all values set to val.
 *
 * @param val Value to initialise the members of the array with.
 * @tparam T Datatype of the array.
 * @tparam size Size of the array.
 */
template <typename T, size_t size> std::array<T, size> make_array(T val) {
  std::array<T, size> ret;
  ret.fill(val);
  return ret;
}

/*!
 * Adds two vectors, taking into account that one vector might be longer than
 * the other one.
 *
 * @param a First array to be added.
 * @param b Second array to be added.
 * @tparam T Datatype of both arrays.
 * @tparam sizea Size of the first array.
 * @tparam sizeb Size of the second array.
 */
template <typename T, size_t sizea, size_t sizeb>
std::array<T, std::max(sizea, sizeb)> add(const std::array<T, sizea> &a,
                                          const std::array<T, sizeb> &b) {
  static constexpr size_t NEW_ARRAY_SIZE = std::max(sizea, sizeb);
  if constexpr (sizeb > sizea) {
    return add(b, a);
  } else {
    std::array<T, NEW_ARRAY_SIZE> ret = a;
    for (size_t i = 0; i < sizeb; i++) {
      ret[i] += b[i];
    }
    return ret;
  }
};

/*!
 * Copies one array into a larger one, filling up the additional members with
 * zeros.
 *
 * @param in Array to be copied.
 * @tparam T Datatype of the arrays.
 * @tparam sizein Size of the input array.
 * @tparam sizeout Size of the output array. Must fulfil sizeout >= sizein.
 */
template <typename T, size_t sizein, size_t sizeout>
std::array<T, sizeout> changearraysize(const std::array<T, sizein> &in) {
  static_assert(sizeout >= sizein,
                "sizeout must be bigger or equal to sizein.");
  if constexpr (sizeout == sizein)
    return in;
  else {
    std::array<T, sizeout> ret = make_array<T, sizeout>(static_cast<T>(0));
    for (size_t i = 0; i < sizein; i++)
      ret[i] = in[i];
    return ret;
  }
}

/*
 * Normally, for every evaluation of a spline, a binary search for the correct
 * interval is necessary. This method is defined in order to integrate every
 * interval separately during the 1D integration, omitting the necessity for the
 * binary search.
 */
template <typename T, size_t ARRAY_SIZE>
T evaluateInterval(const T &x, const std::array<T, ARRAY_SIZE> &coeffs,
                   const T &xm) {
  T result = static_cast<T>(0), xpot = static_cast<T>(1);
  const T dx = x - xm;
  for (const T &c : coeffs) {
    result += c * xpot;
    xpot *= dx;
  }
  return result;
};

}; // end namespace okruz::bspline::internal

#endif // OKRUZ_BSPLINE_MISC_H
