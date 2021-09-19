#ifndef OKRUZ_BSPLINE_OPERATORS_GENERICOPERATORS_H
#define OKRUZ_BSPLINE_OPERATORS_GENERICOPERATORS_H
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
#include <okruz/bspline/exceptions/BSplineException.h>

/*!
 * Operator definitions.
 */
namespace okruz::bspline::operators {
using namespace okruz::bspline::exceptions;

/*!
 * Marker interface for operators. All proper operators must derive from this
 * interface.
 */
class Operator {};

/*!
 * Indicates whether the template parameter is an operator
 * type.
 *
 * @tparam O Template parameter.
 */
template <typename O>
inline constexpr bool is_operator_v = std::is_base_of_v<Operator, O>;

/*!
 * Indicates whether both template parameters are operator
 * types.
 *
 * @tparam O1 First template parameter.
 * @tparam O2 Second template parameter.
 */
template <typename O1, typename O2>
inline constexpr bool are_operators_v = is_operator_v<O1> &&is_operator_v<O2>;

/*!
 * Helper method that applies an operator to a spline based on the
 * transformation of the coefficients on a single interval.
 *
 * @param op The operator to apply to the spline.
 * @param spline The spline to apply the operator to.
 * @tparam T The datatype of the input and output splines.
 * @tparam order The order of the input spline.
 * @tparam O The type of the operator.
 */
template <typename T, size_t order, typename O,
          typename = std::enable_if_t<is_operator_v<O>>>
decltype(auto) transformSpline(O op, const Spline<T, order> &spline) {
  using OutputArray = std::remove_reference_t<std::result_of_t<decltype (
      &O::transform)(const std::array<T, order + 1> &, const T &)>>;

  const auto &oldCoefficients = spline.getCoefficients();

  std::vector<OutputArray> newCoefficients;
  newCoefficients.reserve(oldCoefficients.size());

  for (size_t i = 0; i < oldCoefficients.size(); i++) {
    const T xm = (spline.getSupport().at(i) + spline.getSupport().at(i + 1)) /
                 static_cast<T>(2);
    newCoefficients.push_back(op.transform(oldCoefficients.at(i), xm));
  }

  return Spline(spline.getSupport(), std::move(newCoefficients));
}

/*!
 * Represents the unity operator.
 */
class UnityOperator : public Operator {
 public:
  /*!
   * Applies the operator to a spline.
   *
   * @param spline The spline to apply the operator to.
   * @tparam T The datatype of the spline.
   * @tparam order The order of the spline.
   */
  template <typename T, size_t order>
  Spline<T, order> operator*(const Spline<T, order> &spline) {
    return spline;
  }

  /*!
   * Applies the operator to a set of coefficients (representing a polynomial on
   * one interval).
   *
   * @param input The polynomial coefficients.
   * @param xm The middlepoint of the interval, with respect to wich the
   * polynomial is defined.
   * @tparam T The datatype of the coefficients.
   * @tparam size The size of the array, i. e. the number of coefficients.
   */
  template <typename T, size_t size>
  std::array<T, size> transform(const std::array<T, size> &input,
                                [[maybe_unused]] const T &xm) {
    return input;
  }
};

}  // namespace okruz::bspline::operators
#endif  // OKRUZ_BSPLINE_OPERATORS_GENERICOPERATORS_H
