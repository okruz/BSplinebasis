#ifndef OKRUZ_BSPLINE_INTEGRATION_LINEARFORM_H
#define OKRUZ_BSPLINE_INTEGRATION_LINEARFORM_H
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

#include <okruz/bspline/Spline.h>
#include <okruz/bspline/operators/GenericOperators.h>

namespace okruz::bspline::integration {

/*!
 * Represents a linear form of a BSpline with an operator applied to the spline.
 *
 * @tparam O The type of the operator applied to the spline.
 */
template <typename O,
          std::enable_if_t<operators::is_operator_v<O>, bool> = true>
class LinearForm {
 private:
  /*! Operator applied to the spline.*/
  const O _o;

  /*!
   * Evaluates the linear form on one interval.
   *
   * @param a The coefficients  of the polynomial.
   * @param dxhalf The half width of the interval.
   * @tparam T The datatype of the polynomials.
   * @tparam size The number of coefficients of the polynomial.
   */
  template <typename T, size_t size>
  static T integrate(const std::array<T, size> &a, const T &dxhalf) {
    T result = static_cast<T>(0);
    T power_of_dxhalf = static_cast<T>(2) * dxhalf;
    const T dxhalf_squared = dxhalf * dxhalf;

    for (size_t i = 0; i < size; i += 2) {
      result += power_of_dxhalf * a[i] / static_cast<T>(i + 1);
      power_of_dxhalf *= dxhalf_squared;
    }
    return result;
  }

 public:
  /*!
   * Constructor constructing a LinearForm from an operator.
   *
   * @param o The operator acting on the spline.
   */
  LinearForm(O o) : _o(o){};

  /*!
   * Default constructor constructing a LinearForm.
   */
  LinearForm() : _o(O{}){};

  /*!
   * Evaluates the linear form for a particular spline.
   *
   * @param a The  spline.
   * @tparam T The datatype of the splines.
   * @tparam order The order of the spline.
   */
  template <typename T, size_t order>
  T integrate(const Spline<T, order> &a) const {
    const size_t nintervals = a.getSupport().numberOfIntervals();

    if (nintervals == 0) return static_cast<T>(0);  // no overlap

    T result = static_cast<T>(0);

    for (size_t i = 0; i < nintervals; i++) {
      const T xm =
          (a.getSupport()[i + 1] + a.getSupport()[i]) / static_cast<T>(2);
      const T dxhalf =
          (a.getSupport()[i + 1] - a.getSupport()[i]) / static_cast<T>(2);

      result += integrate(_o.transform(a.getCoefficients()[i], xm), dxhalf);
    }
    return result;
  }
};

/*!
 * Deduction guide for a bilinear form with no operator explicitly defined.
 */
LinearForm()->LinearForm<operators::UnityOperator>;

}  // namespace okruz::bspline::integration
#endif  // OKRUZ_BSPLINE_INTEGRATION_LINEARFORM_H
