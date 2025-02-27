/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#ifndef BSPLINE_OPERATORS_GENERICOPERATORS_H
#define BSPLINE_OPERATORS_GENERICOPERATORS_H

#include <bspline/Spline.h>

/*!
 * @brief Operator definitions.
 */
namespace bspline::operators {

/*!
 * @brief Defines an operator.
 *
 * Indicates whether the template parameter is an operator
 * type.
 */
template <typename O>
concept Operator =
    requires(O o, size_t inputOrder, const std::array<double, 4> &input,
             const support::Grid<double> &grid, size_t intervalIndex) {
      { O::outputOrder(inputOrder) } -> std::same_as<size_t>;
      {
        o.transform(input, grid, intervalIndex)
      } -> std::same_as<std::array<double, O::outputOrder(4 - 1) + 1>>;
    };

/*!
 * @brief Applies operator to spline.
 *
 * Helper method that applies an operator to a spline based on the
 * transformation of the coefficients on a single interval.
 *
 * @param op The operator to apply to the spline.
 * @param spline The spline to apply the operator to.
 * @tparam T The datatype of the input and output splines.
 * @tparam order The order of the input spline.
 * @tparam O The type of the operator.
 * @returns The spline resulting from the application of this operator to the
 * spline.
 */
template <typename T, size_t order, Operator O>
auto transformSpline(const O &op, const Spline<T, order> &spline) {
  constexpr size_t OUTPUT_SIZE = O::outputOrder(order) + 1;

  const auto &oldCoefficients = spline.getCoefficients();

  std::vector<std::array<T, OUTPUT_SIZE>> newCoefficients;
  newCoefficients.reserve(oldCoefficients.size());

  const auto &support = spline.getSupport();

  for (size_t i = 0; i < oldCoefficients.size(); i++) {
    const size_t absIndex = support.absoluteFromRelative(i);
    newCoefficients.push_back(
        op.transform(oldCoefficients[i], support.getGrid(), absIndex));
  }

  return Spline(spline.getSupport(), std::move(newCoefficients));
}

// ################## Generic Operator Definitions #######################
// #######################################################################

// #######################################################################
// ########################## IdentityOperator ##############################

/*!
 * @brief Multiplicative identity operator.
 *
 * Represents the identity operator.
 */
class IdentityOperator final {
 public:
  /*!
   * @brief Order of the resulting Spline.
   *
   * Returns the order of the output Spline for a given input order.
   *
   * @param inputOrder the order of the input Spline.
   * @returns The output Spline-order for a given input input order.
   */
  static constexpr size_t outputOrder(size_t inputOrder) { return inputOrder; }

  /*!
   * @brief Applies operator on one interval.
   *
   * Applies the operator to a set of coefficients (representing a polynomial on
   * one interval).
   *
   * @param input The polynomial coefficients.
   * @param grid The global grid with respect to which the splines are defined.
   * @param intervalIndex The index of the begin of the interval with respect to
   * the global grid.
   * @tparam T The datatype of the coefficients.
   * @tparam size The size of the array, i. e. the number of coefficients.
   * @returns The polyomial coefficients arising from the application of the
   * operator to the input coefficients.
   */
  template <typename T, size_t size>
  std::array<T, size> transform(const std::array<T, size> &input,
                                [[maybe_unused]] const support::Grid<T> &grid,
                                [[maybe_unused]] size_t intervalIndex) const {
    return input;
  }
};

/*!
 * @brief Applies an operator to a spline.
 *
 * @param o The operator.
 * @param s The spline.
 * @tparam O The type of the operator.
 * @tparam T The data type Spline.
 * @tparam order The order of the Spline.
 * @returns The spline resulting from the application of the operator to the
 * spline.
 */
template <Operator O, typename T, size_t order>
auto operator*(const O &o, const Spline<T, order> &s) {
  return transformSpline(o, s);
}
}  // namespace bspline::operators
#endif  // BSPLINE_OPERATORS_GENERICOPERATORS_H
