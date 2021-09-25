#ifndef OKRUZ_BSPLINE_OPERATORS_COMPOUNDOPERATORS_H
#define OKRUZ_BSPLINE_OPERATORS_COMPOUNDOPERATORS_H
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

namespace okruz::bspline::operators {
/*!
 * Represents the product of two operators.
 *
 * @tparam O1 The type of the first (left) operator.
 * @tparam O2 The type of the second (right) operator.
 */
template <typename O1, typename O2,
          std::enable_if_t<are_operators_v<O1, O2>, bool> = true>
class OperatorProduct : public Operator {
 private:
  /*! The first (left) operator.*/
  const O1 _o1;
  /*! The second (right) operator.*/
  const O2 _o2;

 public:
  /*!
   * Creates an OperatorProduct from two operators.
   *
   * @param o1 The first (left) operator.
   * @param o2 The second (right) ooperator.
   */
  OperatorProduct(O1 o1, O2 o2) : _o1(o1), _o2(o2){};

  /*!
   * Returns the order of the output spline for a given input order.
   *
   * @param inputOrder the order of the input spline.
   */
  static constexpr size_t outputOrder(size_t inputOrder) {
    return O1::outputOrder(O2::outputOrder(inputOrder));
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
    return _o1.transform(_o2.transform(input, xm), xm);
  }
};

/*!
 * The multiplication operator for two operators, returning an OperatorProduct.
 *
 * @param o1 The first (left) operator.
 * @param o2 The second (right) operator.
 * @tparam O1 The type of the first (left) operator.
 * @tparam O2 The type of the second (right) operator.
 */
template <typename O1, typename O2,
          std::enable_if_t<are_operators_v<O1, O2>, bool> = true>
OperatorProduct<O1, O2> operator*(const O1 &o1, const O2 &o2) {
  return OperatorProduct(o1, o2);
}

// ######################### OperatorProduct #############################
// #######################################################################

// #######################################################################
// ########################### OperatorSum ###############################

/*! Indicates which operation the OperatorSum is supposed to perform. */
enum class AdditionOperation {
  /*! Perform addition. */
  ADDITION,
  /*! Perform subtraction. */
  SUBTRACTION,
};

/*!
 * Checks whether the operation is a valid oeration for OperatorSum.
 *
 * @tparam operation The operation type to be checked.
 */
template <AdditionOperation operation>
inline constexpr bool is_valid_operation_v = (operation ==
                                              AdditionOperation::ADDITION) ||
                                             (operation ==
                                              AdditionOperation::SUBTRACTION);

/*!
 * Represents the sum or difference of two operators.
 *
 * @tparam O1 The type of the first operator.
 * @tparam O2 The type of the second operator.
 * @tparam operation Indicates whether the operators shall be added or
 * subtracted.
 */
template <
    typename O1, typename O2, AdditionOperation operation,
    std::enable_if_t<are_operators_v<O1, O2> && is_valid_operation_v<operation>,
                     bool> = true>
class OperatorSum : public Operator {
 private:
  /*! The first operator.*/
  const O1 _o1;
  /*! The second operator.*/
  const O2 _o2;

  /*!
   * Adds the smaller array to the larger array in place and returns the larger
   * array.
   *
   * @param a The first array.
   * @param b The second array.
   * @tparam T The datatype of the arrays.
   * @tparam sizea The size of a.
   * @tparam sizeb The size of b.
   */
  template <typename T, size_t sizea, size_t sizeb>
  static std::array<T, std::max(sizea, sizeb)> &add(std::array<T, sizea> &a,
                                                    std::array<T, sizeb> &b) {
    if constexpr (sizeb > sizea) {
      return add(b, a);
    } else {
      for (size_t i = 0; i < sizeb; i++) {
        a[i] += b[i];
      }
      return a;
    }
  }

 public:
  /*!
   * Creates an OperatorSum from two operators.
   *
   * @param o1 The first operator.
   * @param o2 The second ooperator.
   */
  OperatorSum(O1 o1, O2 o2) : _o1(o1), _o2(o2){};

  /*!
   * Returns the order of the output spline for a given input order.
   *
   * @param inputOrder the order of the input spline.
   */
  static constexpr size_t outputOrder(size_t inputOrder) {
    return std::max(O1::outputOrder(inputOrder), O2::outputOrder(inputOrder));
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
    auto a = _o1.transform(input, xm);
    auto b = _o2.transform(input, xm);

    // Negate b if subtraction is requested.
    if constexpr (operation == AdditionOperation::SUBTRACTION) {
      for (T &el : b) {
        el *= static_cast<T>(-1);
      }
    }

    return add(a, b);
  }
};

/*!
 * The addition operator for two operators, returning an OperatorSum.
 *
 * @param o1 The first operator.
 * @param o2 The second operator.
 * @tparam O1 The type of the first operator.
 * @tparam O2 The type of the second operator.
 */
template <typename O1, typename O2,
          std::enable_if_t<are_operators_v<O1, O2>, bool> = true>
auto operator+(const O1 &o1, const O2 &o2) {
  return OperatorSum<O1, O2, AdditionOperation::ADDITION>(o1, o2);
}

/*!
 * The subtraction operator for two operators, returning an OperatorSum.
 *
 * @param o1 The first operator.
 * @param o2 The second operator.
 * @tparam O1 The type of the first operator.
 * @tparam O2 The type of the second operator.
 */
template <typename O1, typename O2,
          std::enable_if_t<are_operators_v<O1, O2>, bool> = true>
auto operator-(const O1 &o1, const O2 &o2) {
  return OperatorSum<O1, O2, AdditionOperation::SUBTRACTION>(o1, o2);
}

}  // namespace okruz::bspline::operators
#endif  // OKRUZ_BSPLINE_OPERATORS_COMPOUNDOPERATORS_H
