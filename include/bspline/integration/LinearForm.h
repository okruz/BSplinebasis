/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#ifndef BSPLINE_INTEGRATION_LINEARFORM_H
#define BSPLINE_INTEGRATION_LINEARFORM_H

#include <bspline/Spline.h>
#include <bspline/operators/GenericOperators.h>

namespace bspline::integration {

/*!
 * Represents the linear form \f[\left\langle
 * a\right\rangle=\int\limits_{-\infty}^{\infty}\mathrm{d}x~\left[\hat{O}\,\,
 * a(x)\right]\f] with the operator \f$\hat{O}\f$ applied to the spline.
 *
 * @tparam O The type of the operator applied to the spline.
 */
template <typename O,
          std::enable_if_t<operators::is_operator_v<O>, bool> = true>
class LinearForm final {
 private:
  /*! Operator applied to the spline.*/
  O _o;

  /*!
   * Evaluates the linear form on one interval.
   *
   * @param a The coefficients  of the polynomial.
   * @param dxhalf The half width of the interval.
   * @tparam T The datatype of the polynomials.
   * @tparam size The number of coefficients of the polynomial.
   * @returns The value of the linear form on the one interval.
   */
  template <typename T, size_t size>
  static T evaluateInterval(const std::array<T, size> &a, const T &dxhalf) {
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
  explicit LinearForm(O o) : _o(std::move(o)){};

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
   * @returns The value of the linear form for the given spline.
   */
  template <typename T, size_t order>
  T evaluate(const Spline<T, order> &a) const {
    const size_t nintervals = a.getSupport().numberOfIntervals();

    T result = static_cast<T>(0);

    for (size_t i = 0; i < nintervals; i++) {
      const size_t absIndex = a.getSupport().absoluteFromRelative(i);
      const T dxhalf =
          (a.getSupport()[i + 1] - a.getSupport()[i]) / static_cast<T>(2);

      result +=
          evaluateInterval(_o.transform(a.getCoefficients()[i],
                                        a.getSupport().getGrid(), absIndex),
                           dxhalf);
    }
    return result;
  }
};

/*!
 * Deduction guide for a linear form which represents the integral over the
 * spline \f[\left\langle
 * a\right\rangle=\int\limits_{-\infty}^{\infty}\mathrm{d}x~a(x).\f]
 */
LinearForm()->LinearForm<operators::IdentityOperator>;

}  // namespace bspline::integration
#endif  // BSPLINE_INTEGRATION_LINEARFORM_H
