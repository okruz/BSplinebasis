#ifndef BSPLINE_OPERATORS_SCALAROPERATORS_H
#define BSPLINE_OPERATORS_SCALAROPERATORS_H
/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <bspline/Spline.h>
#include <bspline/operators/CompoundOperators.h>

/*!
 * Operator definitions.
 */
namespace bspline::operators {

/*!
 * Indicates whether the argument types may be arguments of a scalar
 * multiplication (i.e. a scalar multiplied with an operator).
 *
 * @tparam T Scalar type.
 * @tparam O Operator type.
 */
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
  ScalarMultiplication(S s, O o) : _s(std::move(s)), _o(std::move(o)){};

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
ScalarMultiplication<S, O> operator*(const S &s, O &&o) {
  return ScalarMultiplication(s, std::forward<O>(o));
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
ScalarMultiplication<S, O> operator*(O &&o, const S &s) {
  return ScalarMultiplication(s, std::forward<O>(o));
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
ScalarMultiplication<S, O> operator/(O &&o, const S &s) {
  return ScalarMultiplication(static_cast<S>(1) / s, std::forward<O>(o));
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
auto operator+(O &&o, const S &s) {
  return std::forward<O>(o) + ScalarMultiplication{s};
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
auto operator+(const S &s, O &&o) {
  return ScalarMultiplication{s} + std::forward<O>(o);
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
auto operator-(O &&o, const S &s) {
  return std::forward<O>(o) - ScalarMultiplication{s};
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
auto operator-(const S &s, O &&o) {
  return ScalarMultiplication{s} - std::forward<O>(o);
}

/*!
 * The unitary minus operator for an operator. Returns an
 * ScalarMultiplication<int, O>.
 *
 * @param o The operator to be negated.
 * @tparam O The type of the operator to be negated.
 */
template <typename O, std::enable_if_t<is_operator_v<O>, bool> = true>
ScalarMultiplication<int, O> operator-(O &&o) {
  return ScalarMultiplication<int, O>(-1, std::forward<O>(o));
}
}  // namespace bspline::operators
#endif  // BSPLINE_OPERATORS_SCALAROPERATORS_H
