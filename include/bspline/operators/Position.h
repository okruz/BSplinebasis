#ifndef BSPLINE_OPERATORS_POSITION_H
#define BSPLINE_OPERATORS_POSITION_H
/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <bspline/internal/misc.h>
#include <bspline/operators/GenericOperators.h>

/*!
 * Operator definitions.
 */
namespace bspline::operators {
namespace internal = bspline::internal;

/*!
 * Represents a power of the position operator \f$x^n\f$.
 *
 * @tparam n Order of the power.
 */
template <size_t n>
class Position : public Operator {
 private:
  /*!
   * Expands \f$(x + x_m)^n\f$.
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
 * Alias for the power \f$x^n\f$ of the position operator.
 *
 * @tparam n Order of the power of the position operator.
 */
template <size_t n>
using X = Position<n>;

}  // namespace bspline::operators
#endif  // BSPLINE_OPERATORS_POSITION_H
