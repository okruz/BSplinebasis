#ifndef OKRUZ_BSPLINE_OPERATORS_POSITION_H
#define OKRUZ_BSPLINE_OPERATORS_POSITION_H
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

#include <okruz/bspline/internal/misc.h>
#include <okruz/bspline/operators/GenericOperators.h>

/*!
 * Operator definitions.
 */
namespace okruz::bspline::operators {
namespace internal = okruz::bspline::internal;

/*!
 * Represents a power of the position operator x^n.
 *
 * @tparam n Order of the power.
 */
template <size_t n>
class Position : public Operator {
 private:
  /*!
   * Expands (x + xm)^n
   *
   * @param xm The middle point of the interval.
   * @tparam T The datatype of xm and the returned coefficients.
   */
  template <typename T>
  static std::array<T, n + 1> expandPower(const T &xm) {
    std::array<T, n + 1> retVal;
    T power_of_xm = static_cast<T>(1);
    for (size_t i = 0; i < n + 1; i++) {
      retVal[n - i] = internal::binomialCoefficient<T>(n, i) * power_of_xm;
      power_of_xm *= xm;
    }
    return retVal;
  }

 public:
  /*!
   * Returns the order of the output spline for a given input order.
   *
   * @param inputOrder the order of the input spline.
   */
  static constexpr size_t outputOrder(size_t inputOrder) {
    return inputOrder + n;
  }

  /*!
   * Applies the operator to a spline.
   *
   * @param spline The spline to apply the operator to.
   * @tparam T The datatype of the splines.
   * @tparam order The order of the input spline.
   */
  template <typename T, size_t order>
  Spline<T, outputOrder(order)> operator*(
      const Spline<T, order> &spline) const {
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
      const std::array<T, size> &input, const T &xm) const {
    constexpr size_t OUTPUT_SIZE = size + n;

    const std::array<T, n + 1> expanded = expandPower<T>(xm);

    std::array<T, OUTPUT_SIZE> retVal;
    retVal.fill(static_cast<T>(0));

    for (size_t i = 0; i < input.size(); i++) {
      for (size_t j = 0; j < expanded.size(); j++) {
        retVal[i + j] += expanded[j] * input[i];
      }
    }
    return retVal;
  }
};

/*!
 * Alias for the position opertor.
 *
 * @tparam n Order of the power of the position operator.
 */
template <size_t n>
using X = Position<n>;

}  // namespace okruz::bspline::operators
#endif  // OKRUZ_BSPLINE_OPERATORS_POSITION_H
