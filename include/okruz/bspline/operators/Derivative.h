#ifndef OKRUZ_BSPLINE_OPERATORS_DERIVATIVE_H
#define OKRUZ_BSPLINE_OPERATORS_DERIVATIVE_H
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

#include <okruz/bspline/operators/GenericOperators.h>

/*!
 * Operator definitions.
 */
namespace okruz::bspline::operators {

/*!
 * Represents a derivative operator of order n, i.e. d^n/dx^n.
 *
 * @tparam n Order of the derivative.
 */
template <size_t n = 1>
class Derivative : public Operator {
 private:
  /*!
   * returns the faculty ratio counter! / denominator!.
   *
   * @param counter The faculty of this value is the counter.
   * @param denominator The faculty of this value is the denominator.
   * @tparam T The datatype of the return type.
   */
  template <typename T>
  static T facultyRatio(size_t counter, size_t denominator) {
    if (denominator > counter) {
      return static_cast<T>(1) / facultyRatio<T>(denominator, counter);
    } else {
      T retVal = static_cast<T>(1);
      for (size_t i = denominator + 1; i <= counter; i++) {
        retVal *= static_cast<T>(i);
      }
      return retVal;
    }
  }

 public:
  /*!
   * Returns the order of the output spline for a given input order.
   *
   * @param inputOrder the order of the input spline.
   */
  static constexpr size_t outputOrder(size_t inputOrder) {
    return std::max(n, inputOrder) - n;
  }

  /*!
   * Applies the operator to a spline.
   *
   * @param spline The spline to apply the operator to.
   * @tparam T The datatype of the splines.
   * @tparam order The order of the input spline.
   */
  template <typename T, size_t order>
  Spline<T, outputOrder(order)> operator*(const Spline<T, order> &spline) {
    return transformSpline(*this, spline);
  }

  /*!
   * Applies the operator to a set of coefficients (representing a polynomial on
   * one interval).
   *
   * @param input The polynomial coefficients.
   * @param xm The middlepoint of the interval, with respect to which the
   * polynomial is defined.
   * @tparam T The datatype of the coefficients.
   * @tparam size The size of the input array, i. e. the number of coefficients.
   */
  template <typename T, size_t size>
  std::array<T, outputOrder(size - 1) + 1> transform(
      const std::array<T, size> &input, [[maybe_unused]] const T &xm) {
    static_assert(size >= 1, "Arrays of size zero not supported.");
    // The order of the input spline.
    constexpr size_t SPLINE_ORDER = size - 1;
    // The size of the output array.
    constexpr size_t OUTPUT_SIZE = outputOrder(SPLINE_ORDER) + 1;

    if constexpr (n > SPLINE_ORDER) {
      return {static_cast<T>(0)};
    } else {
      std::array<T, OUTPUT_SIZE> retVal;
      for (size_t i = 0; i < OUTPUT_SIZE; i++) {
        retVal[i] = facultyRatio<T>(i + n, i) * input[i + n];
      }
      return retVal;
    }
  }
};

/*!
 * Alias for the derivative opertor.
 *
 * @tparam n Order of the derivative.
 */
template <size_t n = 1>
using Dx = Derivative<n>;

}  // namespace okruz::bspline::operators
#endif  // OKRUZ_BSPLINE_OPERATORS_DERIVATIVE_H
