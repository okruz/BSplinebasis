#ifndef OKRUZ_BSPLINE_OPERATORS_COMPOUNDOPERATOR_H
#define OKRUZ_BSPLINE_OPERATORS_COMPOUNDOPERATOR_H
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

#include <okruz/bspline/operators/GenericOperator.h>

namespace okruz::bspline::operators {

/*!
 * Convenience struct, to determine if both template parameters are operator
 * types.
 *
 * @tparam O1 First template parameter.
 * @tparam O2 Second template parameter.
 */
template <typename O1, typename O2>
struct are_operators {
  /*! Indicates whether both template parameter types are operator types.*/
  static constexpr bool value = std::is_base_of<Operator, O1>::value &&
                                std::is_base_of<Operator, O2>::value;
};

/*!
 * A compound operator represents the product of two operators.
 *
 * @tparam O1 The type of the first (left) operator.
 * @tparam O2 The type of the second (right) operator.
 */
template <typename O1, typename O2,
          typename = std::enable_if_t<are_operators<O1, O2>::value>>
class CompoundOperator : public Operator {
 private:
  /*! The first (left) operator.*/
  O1 _o1;
  /*! The second (right) operator.*/
  O2 _o2;

 public:
  /*!
   * Creates a compound operator as a product of two operators.
   *
   * @param o1 The first (left) operator.
   * @param o2 The second (right) ooperator.
   */
  CompoundOperator(O1 o1, O2 o2) : m_o1(o1), m_o2(o2){};

  /*!
   * Applies the operator to a spline.
   *
   * @param spline The spline to apply the operator to.
   * @tparam T The datatype of the splines.
   * @tparam order The order of the input spline.
   */
  template <typename T, size_t order>
  decltype(auto) operator*(const Spline<T, order> &spline) {
    return transformSpline(*this, spline);
  }

  /*!
   * Applies the operator to a set of coefficients (representing a polynomial on
   * one interval).
   *
   * @param input The polynomial coefficients.
   * @param xm The middlepoint of the interval, with respect to wich the
   * polynomial is defined.
   * @tparam T The datatype of the coefficients.
   * @tparam size The size of the input array, i. e. the number of coefficients.
   */
  template <typename T, size_t size>
  decltype(auto) transform(const std::array<T, size> &input, const T &xm) {
    return _o1.transform(_o2.transform(input), xm), xm);
  }
};

/*!
 * The multiplication operator for two operators, returning a CompoundOperator.
 *
 * @tparam O1 The type of the first (left) operator.
 * @tparam O2 The type of the second (right) operator.
 */
template <typename O1, typename O2,
          typename = std::enable_if_t<are_operators<O1, O2>::value>>
CompoundOperator<O1, O2> operator*(const O1 &o1, const O2 &o2) {
  return CompoundOperator(o1, o2);
};

}  // namespace okruz::bspline::operators
#endif  // OKRUZ_BSPLINE_OPERATORS_COMPOUNDOPERATOR_H
