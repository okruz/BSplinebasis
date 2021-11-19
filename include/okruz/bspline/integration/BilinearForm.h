#ifndef OKRUZ_BSPLINE_INTEGRATION_BILINEARFORM_H
#define OKRUZ_BSPLINE_INTEGRATION_BILINEARFORM_H
/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <okruz/bspline/Spline.h>
#include <okruz/bspline/operators/GenericOperators.h>

namespace okruz::bspline::integration {

/*!
 * Represents a bilinear form of two BSplines with two different operators
 * applied to the splines.
 *
 * @tparam O1 The type of the operator applied to the first spline.
 * @tparam O2 The type of the operator applied to the second spline.
 */
template <typename O1, typename O2,
          std::enable_if_t<operators::are_operators_v<O1, O2>, bool> = true>
class BilinearForm {
 private:
  /*! Operator applied to the first spline.*/
  const O1 _o1;
  /*! Operator applied to the second spline.*/
  const O2 _o2;

  /*!
   * Evaluates the bilinear form on one interval.
   *
   * @param a The coefficients  of the first polynomial.
   * @param b The coefficients of the second polynomial.
   * @param dxhalf The half width of the interval.
   * @tparam T The datatype of the polynomials.
   * @tparam sizea The number of coefficients of the first polynomial.
   * @tparam sizeb The number of coefficients of the second polynomial.
   */
  template <typename T, size_t sizea, size_t sizeb>
  static T integrate(const std::array<T, sizea> &a,
                     const std::array<T, sizeb> &b, const T &dxhalf) {
    std::array<T, (sizea + sizeb) / 2> coefficients;
    coefficients.fill(static_cast<T>(0));
    for (size_t i = 0; i < sizea; i++) {
      for (size_t j = i % 2; j < sizeb; j += 2) {
        coefficients[(i + j) / 2] += a[i] * b[j];
      }
    }

    T result = static_cast<T>(0);
    T power_of_dxhalf = static_cast<T>(2) * dxhalf;
    const T dxhalf_squared = dxhalf * dxhalf;

    for (size_t i = 0; i < coefficients.size(); i++) {
      const size_t order = i * 2;
      result += power_of_dxhalf * coefficients[i] / static_cast<T>(order + 1);
      power_of_dxhalf *= dxhalf_squared;
    }
    return result;
  }

 public:
  /*!
   * Constructor constructing a BilinearForm from two operators.
   *
   * @param o1 The operator acting on the first (left) spline.
   * @param o2 The operator acting on the second (right) spline.
   */
  BilinearForm(O1 o1, O2 o2) : _o1(std::move(o1)), _o2(std::move(o2)){};

  /*!
   * Constructor constructing a BilinearForm from one operators.
   *
   * @param o2 The operator acting on the second (right) spline.
   */
  explicit BilinearForm(O2 o2) : _o1(O1{}), _o2(std::move(o2)){};

  /*!
   * Default constructor constructing a BilinearForm.
   */
  BilinearForm() : _o1(O1{}), _o2(O2{}){};

  /*!
   * Evaluates the bilinear form for two particular splines.
   *
   * @param a The first (left) spline.
   * @param b The second (right) spline.
   * @tparam T The datatype of the splines.
   * @tparam ordera The order of the first (left) spline.
   * @tparam orderb The order of the second (right) spline.
   */
  template <typename T, size_t ordera, size_t orderb>
  T integrate(const Spline<T, ordera> &a, const Spline<T, orderb> &b) const {
    // Will also check whether the two grids are equivalent.
    const support::Support integrandSupport =
        a.getSupport().calcIntersection(b.getSupport());
    const size_t nintervals = integrandSupport.numberOfIntervals();

    T result = static_cast<T>(0);

    for (size_t interv = 0; interv < nintervals; interv++) {
      const auto absIndex = integrandSupport.absoluteFromRelative(interv);

      const auto aIndex =
          a.getSupport().intervalIndexFromAbsolute(absIndex).value();
      const auto bIndex =
          b.getSupport().intervalIndexFromAbsolute(absIndex).value();

      const T xm = (a.getSupport()[aIndex + 1] + a.getSupport()[aIndex]) /
                   static_cast<T>(2);
      const T dxhalf = (a.getSupport()[aIndex + 1] - a.getSupport()[aIndex]) /
                       static_cast<T>(2);

      result +=
          integrate(_o1.transform(a.getCoefficients()[aIndex], xm),
                    _o2.transform(b.getCoefficients()[bIndex], xm), dxhalf);
    }
    return result;
  }
};

/*!
 * Deduction guide for a bilinear form with only one operator to be applied to
 * the second (right) spline.
 *
 * @tparam O2 Type of the operator applied to the second spline.
 */
template <typename O2>
BilinearForm(O2 o2) -> BilinearForm<operators::UnityOperator, O2>;

/*!
 * Deduction guide for a bilinear form with no operator explicitly defined.
 */
BilinearForm()
    ->BilinearForm<operators::UnityOperator, operators::UnityOperator>;

/*!
 * Short hand for a scalar product.
 */
using ScalarProduct =
    BilinearForm<operators::UnityOperator, operators::UnityOperator>;

}  // namespace okruz::bspline::integration
#endif  // OKRUZ_BSPLINE_INTEGRATION_BILINEARFORM_H
