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
#include <okruz/bspline/operators/CompoundOperators.h>

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
 * @tparam S type of the scalar.
 * @tparam O type of the operator to be multiplied with the scalar.
 */
template <typename S, typename O,
          std::enable_if_t<are_scalar_multiplication_types_v<S, O>, bool> =
              true>
class ScalarMultiplication : public Operator {
 private:
  /*! The scalar to multiply the operator with. */
  const S _s;
  /*! The operator to be multiplied. */
  const O _o;

 public:
  /*!
   * Constructor constructing an ScalarMultiplication from a scalar and an
   * operator.
   *
   * @param s The scalar to be multiplied.
   * @param o The operator to be multiplied.
   */
  ScalarMultiplication(S s, O o) : _s(s), _o(o){};

  /*!
   * Constructor constructing an ScalarMultiplication from a scalar alone.
   *
   * @param s The scalar to be multiplied.
   */
  ScalarMultiplication(S s) : _s(s), _o(UnityOperator{}){};

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
  auto operator*(const Spline<S, order> &spline) const {
    return transformSpline(*this, spline);
  }

  /*!
   * Applies the operator to a set of coefficients (representing a polynomial on
   * one interval).
   *
   * @param input The polynomial coefficients.
   * @param xm The middlepoint of the interval, with respect to which the
   * polynomial is defined.
   * @tparam T The datatype of the spline.
   * @tparam size The size of the array, i. e. the number of coefficients.
   */
  template <typename T, size_t size>
  auto transform(const std::array<T, size> &input, const T &xm) const {
    auto a = _o.transform(input, xm);

    // Multiply a.
    for (T &el : a) {
      el *= static_cast<T>(_s);
    }
    return a;
  }
};

/*!
 * Deduction guide for a ScalarMultiplication constructed from a scalar alone.
 *
 * @tparam S The type of the scalar.
 */
template <typename S>
ScalarMultiplication(S s) -> ScalarMultiplication<S, UnityOperator>;

/*!
 * The scalar multiplication operator for an operator.
 *
 * @param s The scalar to be multiplied.
 * @param o The operator to be multiplied.
 * @tparam S The type of the scalar.
 * @tparam O The type of the operator to be multiplied.
 */
template <
    typename S, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<S, O>, bool> = true>
ScalarMultiplication<S, O> operator*(const S &s, const O &o) {
  return ScalarMultiplication(s, o);
}

/*!
 * The scalar multiplication operator for an operator.
 *
 * @param o The operator to be multiplied.
 * @param s The scalar to be multiplied.
 * @tparam S The type of the scalar.
 * @tparam O The type of the operator to be multiplied.
 */
template <
    typename S, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<S, O>, bool> = true>
ScalarMultiplication<S, O> operator*(const O &o, const S &s) {
  return ScalarMultiplication(s, o);
}

/*!
 * The scalar division operator for an operator.
 *
 * @param o The operator to be divided.
 * @param s The divisor.
 * @tparam S The type of the scalar.
 * @tparam O The type of the operator to be divided.
 */
template <
    typename S, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<S, O>, bool> = true>
ScalarMultiplication<S, O> operator/(const O &o, const S &s) {
  return ScalarMultiplication(static_cast<S>(1) / s, o);
}

/*!
 * The scalar addition operator for an operator.
 *
 * @param o The operator to be added.
 * @param s The scalar to be added.
 * @tparam S The type of the scalar.
 * @tparam O The type of the operator.
 */
template <
    typename S, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<S, O>, bool> = true>
auto operator+(const O &o, const S &s) {
  return o + ScalarMultiplication{s};
}

/*!
 * The scalar addition operator for an operator.
 *
 * @param s The scalar to be added.
 * @param o The operator to be added.
 * @tparam S The type of the scalar.
 * @tparam O The type of the operator.
 */
template <
    typename S, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<S, O>, bool> = true>
auto operator+(const S &s, const O &o) {
  return ScalarMultiplication{s} + o;
}

/*!
 * The scalar subtraction operator for an operator.
 *
 * @param o The operator.
 * @param s The scalar to be subtracted.
 * @tparam S The type of the scalar.
 * @tparam O The type of the operator.
 */
template <
    typename S, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<S, O>, bool> = true>
auto operator-(const O &o, const S &s) {
  return o - ScalarMultiplication{s};
}

/*!
 * The scalar subtraction operator for an operator.
 *
 * @param s The scalar.
 * @param o The operator to be subtracted.
 * @tparam S The type of the scalar.
 * @tparam O The type of the operator.
 */
template <
    typename S, typename O,
    std::enable_if_t<are_scalar_multiplication_types_v<S, O>, bool> = true>
auto operator-(const S &s, const O &o) {
  return ScalarMultiplication{s} - o;
}

/*!
 * The unitary minus operator for an operator. Returns an
 * ScalarMultiplication<int, O>.
 *
 * @param o The operator to be negated.
 * @tparam O The type of the operator to be negated.
 */
template <typename O, std::enable_if_t<is_operator_v<O>, bool> = true>
ScalarMultiplication<int, O> operator-(const O &o) {
  return ScalarMultiplication<int, O>(-1, o);
}
}  // namespace okruz::bspline::operators
#endif  // OKRUZ_BSPLINE_OPERATORS_SCALAROPERATORS_H
