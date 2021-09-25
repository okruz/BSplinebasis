#ifndef OKRUZ_BSPLINE_OPERATORS_SCALAROPERATORS_H
#define OKRUZ_BSPLINE_OPERATORS_SCALAROPERATORS_H
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

/*!
 * Operator definitions.
 */
namespace okruz::bspline::operators {

template <typename T, typename O>
inline constexpr bool are_scalar_multiplication_types_v =
    !is_operator_v<T> && is_operator_v<O>;

/*!
 * Represents the multiplication of an operator with a scalar.
 *
 * @tparam T type of the scalar.
 * @tparam O type of the operator to be multiplied with the scalar.
 */
template <typename T, typename O,
          std::enable_if_t<are_scalar_multiplication_types_v<T, O>, bool> =
              true>
class ScalarMultiplication : public Operator {
 private:
  /*! The scalar to multiply the operator with. */
  T _t;
  /*! The operator to be multiplied. */
  O _o;

 public:
  /*!
   * Constructor constructing an ScalarMultiplication from a scalar and an
   * operator.
   *
   * @param t The scalar to be multiplied.
   * @param o The operator to be multiplied.
   */
  ScalarMultiplication(T t, O o) : _t(t), _o(o){};

  /*!
   * Returns the order of the output spline for a given input order.
   *
   * @param inputOrder the order of the input spline.
   */
  static constexpr size_t outputOrder(size_t inputOrder) {
    return O::outputOrder(inputOrder);
  }

  /*!
   * Applies the operator to a spline.
   *
   * @param spline The spline to apply the operator to.
   * @tparam order The order of the spline.
   */
  template <size_t order>
  auto operator*(const Spline<T, order> &spline) {
    return transformSpline(*this, spline);
  }

  /*!
   * Applies the operator to a set of coefficients (representing a polynomial on
   * one interval).
   *
   * @param input The polynomial coefficients.
   * @param xm The middlepoint of the interval, with respect to which the
   * polynomial is defined.
   * @tparam size The size of the array, i. e. the number of coefficients.
   */
  template <size_t size>
  auto transform(const std::array<T, size> &input, const T &xm) {
    auto a = _o.transform(input, xm);

    // Negate a.
    for (T &el : a) {
      el *= _t;
    }
    return a;
  }
};

/*!
 * The scalar multiplication operator for an operator.
 *
 * @param t The scalar to be multiplied.
 * @param o The operator to be multiplied.
 * @tparam T The type of the scalar.
 * @tparam O The type of the operator to be multiplied.
 */
template <
    typename T, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<T, O>, bool> = true>
ScalarMultiplication<T, O> operator*(const T &t, const O &o) {
  return ScalarMultiplication(t, o);
}

/*!
 * The scalar multiplication operator for an operator.
 *
 * @param o The operator to be multiplied.
 * @param t The scalar to be multiplied.
 * @tparam T The type of the scalar.
 * @tparam O The type of the operator to be multiplied.
 */
template <
    typename T, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<T, O>, bool> = true>
ScalarMultiplication<T, O> operator*(const O &o, const T &t) {
  return ScalarMultiplication(t, o);
}

/*!
 * The scalar division operator for an operator.
 *
 * @param o The operator to be divided.
 * @param t The divisor.
 * @tparam T The type of the scalar.
 * @tparam O The type of the operator to be divided.
 */
template <
    typename T, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<T, O>, bool> = true>
ScalarMultiplication<T, O> operator/(const O &o, const T &t) {
  return ScalarMultiplication(static_cast<T>(1) / t, o);
}

// ########################## UnityOperator ##############################
// #######################################################################

// #######################################################################
// ######################### OperatorNegation ############################

/*!
 * Represents the negation of an operator.
 *
 * @tparam O type of the operator to be negated.
 */
template <typename O, std::enable_if_t<is_operator_v<O>, bool> = true>
class OperatorNegation : public Operator {
 private:
  /*! The operator to be negated. */
  O _o;

 public:
  /*!
   * Constructor constructing an OperatorNegation from an operator.
   *
   * @param o The operator to be negated.
   */
  OperatorNegation(O o) : _o(o){};

  /*!
   * Returns the order of the output spline for a given input order.
   *
   * @param inputOrder the order of the input spline.
   */
  static constexpr size_t outputOrder(size_t inputOrder) {
    return O::outputOrder(inputOrder);
  }

  /*!
   * Applies the operator to a spline.
   *
   * @param spline The spline to apply the operator to.
   * @tparam T The datatype of the spline.
   * @tparam order The order of the spline.
   */
  template <typename T, size_t order>
  auto operator*(const Spline<T, order> &spline) {
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
   * @tparam size The size of the array, i. e. the number of coefficients.
   */
  template <typename T, size_t size>
  auto transform(const std::array<T, size> &input, const T &xm) {
    auto a = _o.transform(input, xm);

    // Negate a.
    for (T &el : a) {
      el *= static_cast<T>(-1);
    }
    return a;
  }
};

/*!
 * The unitary minus operator for an operator. Returns an OperatorNegation.
 *
 * @param o The operator to be negated.
 * @tparam O The type of the operator to be negated.
 */
template <typename O, std::enable_if_t<is_operator_v<O>, bool> = true>
OperatorNegation<O> operator-(const O &o) {
  return OperatorNegation(o);
}
}  // namespace okruz::bspline::operators
#endif  // OKRUZ_BSPLINE_OPERATORS_SCALAROPERATORS_H
