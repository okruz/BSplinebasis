#ifndef BSPLINE_INTEGRATION_BILINEARFORM_H
#define BSPLINE_INTEGRATION_BILINEARFORM_H
/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <bspline/Spline.h>
#include <bspline/operators/GenericOperators.h>

/*!
 * Nampespace containing the integration routines. Analytical integration is
 * represented by the linear and bilinear forms. There is also code for
 * numerical integration using the boost fixed poind Gauss-Legendre scheme.
 */
namespace bspline::integration {

/*!
 * Represents the bilinear form \f[\left\langle a,\, b\right\rangle =
 * \left\langle \hat{O}_1\,a\,\middle|\,\hat{O}_2\,b\right\rangle =
 * \int\limits_{-\infty}^{\infty} \mathrm{d}x~\left[\hat{O}_1\,a(x)\right]
 * \,\,\left[\hat{O}_2\,b(x)\right] \f] with the operators
 * \f$\hat{O}_1,\,\hat{O}_2\f$ applied to the two splines.
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
  static T evaluateInterval(const std::array<T, sizea> &a,
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
   * Constructor constructing a BilinearForm from the two operatos.
   * @param o1 The operator acting on the first (left) spline.
   * @param o2 The operator acting on the second (right) spline.
   */
  BilinearForm(O1 o1, O2 o2) : _o1(std::move(o1)), _o2(std::move(o2)){};

  /*!
   * Constructor constructing a BilinearForm from only the \f$\hat{O}_2\f$
   * operator. The \f$\hat{O}_1\f$-operator is default constructed.
   *
   * @param o2 The operator acting on the second (right) spline.
   */
  explicit BilinearForm(O2 o2) : _o1(O1{}), _o2(std::move(o2)){};

  /*!
   * Default constructor which default constructs both operators.
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
  T evaluate(const Spline<T, ordera> &a, const Spline<T, orderb> &b) const {
    // Will also check whether the two grids are equivalent.
    const support::Support integrandSupport =
        a.getSupport().calcIntersection(b.getSupport());
    const size_t nintervals = integrandSupport.numberOfIntervals();

    const auto &grid = integrandSupport.getGrid();

    T result = static_cast<T>(0);

    for (size_t interv = 0; interv < nintervals; interv++) {
      const auto absIndex = integrandSupport.absoluteFromRelative(interv);

      const auto aIndex =
          a.getSupport().intervalIndexFromAbsolute(absIndex).value();
      const auto bIndex =
          b.getSupport().intervalIndexFromAbsolute(absIndex).value();

      const T dxhalf = (a.getSupport()[aIndex + 1] - a.getSupport()[aIndex]) /
                       static_cast<T>(2);

      result += evaluateInterval(
          _o1.transform(a.getCoefficients()[aIndex], grid, absIndex),
          _o2.transform(b.getCoefficients()[bIndex], grid, absIndex), dxhalf);
    }
    return result;
  }
};

/*!
 * Deduction guide for a bilinear form which corresponds to the matrix element
 * of the operator \f$\hat{O}_2\f$. \f[\left\langle a,\, b\right\rangle  =
 * \left\langle a\,\middle|\,\hat{O}_2\,\middle|\,b\right\rangle =
 * \int\limits_{-\infty}^{\infty} \mathrm{d}x~a(x) \,\,\hat{O}_2\,\,b(x)\ \f].
 *
 * @tparam O2 Type of the operator applied to the second spline.
 */
template <typename O2>
BilinearForm(O2 o2) -> BilinearForm<operators::IdentityOperator, O2>;

/*!
 * Deduction guide for a bilinear form which corresponds to the scalar product.
 * \f[\left\langle a,\, b\right\rangle = \left\langle
 * a\,\middle|\,b\right\rangle = \int\limits_{-\infty}^{\infty} \mathrm{d}x~a(x)
 * \,\,b(x)\ \f]
 */
BilinearForm()
    ->BilinearForm<operators::IdentityOperator, operators::IdentityOperator>;

/*!
 * Short hand for a scalar product \f[\left\langle a,\, b\right\rangle =
 * \left\langle a\,\middle|\,b\right\rangle =  \int\limits_{-\infty}^{\infty}
 * \mathrm{d}x~a(x) \,\,b(x).\ \f]
 */
using ScalarProduct =
    BilinearForm<operators::IdentityOperator, operators::IdentityOperator>;

}  // namespace bspline::integration
#endif  // BSPLINE_INTEGRATION_BILINEARFORM_H
