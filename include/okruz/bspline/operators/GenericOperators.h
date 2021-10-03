#ifndef OKRUZ_BSPLINE_OPERATORS_GENERICOPERATORS_H
#define OKRUZ_BSPLINE_OPERATORS_GENERICOPERATORS_H
/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <okruz/bspline/Spline.h>

/*!
 * Operator definitions.
 */
namespace okruz::bspline::operators {

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
          std::enable_if_t<is_operator_v<O>, bool> = true>
decltype(auto) transformSpline(const O &op, const Spline<T, order> &spline) {
  constexpr size_t OUTPUT_SIZE = O::outputOrder(order) + 1;

  const auto &oldCoefficients = spline.getCoefficients();

  std::vector<std::array<T, OUTPUT_SIZE>> newCoefficients;
  newCoefficients.reserve(oldCoefficients.size());

  for (size_t i = 0; i < oldCoefficients.size(); i++) {
    const T xm = (spline.getSupport().at(i) + spline.getSupport().at(i + 1)) /
                 static_cast<T>(2);
    newCoefficients.push_back(op.transform(oldCoefficients.at(i), xm));
  }

  return Spline(spline.getSupport(), std::move(newCoefficients));
}

// ################## Generic Operator Definitions #######################
// #######################################################################

// #######################################################################
// ########################## UnityOperator ##############################

/*!
 * Represents the unity operator.
 */
class UnityOperator : public Operator {
 public:
  /*!
   * Returns the order of the output spline for a given input order.
   *
   * @param inputOrder the order of the input spline.
   */
  static constexpr size_t outputOrder(size_t inputOrder) { return inputOrder; }

  /*!
   * Applies the operator to a spline.
   *
   * @param spline The spline to apply the operator to.
   * @tparam T The datatype of the spline.
   * @tparam order The order of the spline.
   */
  template <typename T, size_t order>
  Spline<T, order> operator*(const Spline<T, order> &spline) const {
    return spline;
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
  std::array<T, size> transform(const std::array<T, size> &input,
                                [[maybe_unused]] const T &xm) const {
    return input;
  }
};
}  // namespace okruz::bspline::operators
#endif  // OKRUZ_BSPLINE_OPERATORS_GENERICOPERATORS_H
