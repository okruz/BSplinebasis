#ifndef BSPLINE_OPERATORS_SPLINEOPERATOR_H
#define BSPLINE_OPERATORS_SPLINEOPERATOR_H
/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <bspline/Spline.h>
#include <bspline/exceptions/BSplineException.h>
#include <bspline/internal/misc.h>
#include <bspline/operators/GenericOperators.h>

namespace bspline::operators {

/*!
 * Implements the operator representation of a Spline.
 *
 * @tparam T The data type of the spline.
 * @tparam order The order of the spline.
 */
template <typename T, size_t order>
class SplineOperator : public Operator {
 private:
  /*! The spline which this operator represents. */
  const Spline<T, order> _s;

 public:
  /*!
   * Constructor constructing a SplineOperator from a Spline.
   *
   * @param s The spline.
   */
  SplineOperator(Spline<T, order> s) : _s(std::move(s)){};

  /*!
   * Returns the order of the output spline for a given input order.
   *
   * @param inputOrder the order of the input spline.
   */
  static constexpr size_t outputOrder(size_t inputOrder) {
    return inputOrder + order;
  }

  /*!
   * Applies the operator to a spline.
   *
   * @param spline The spline to apply the operator to.
   * @tparam orderSpline The order of the spline.
   */
  template <size_t orderSpline>
  auto operator*(const Spline<T, orderSpline> &spline) const {
    return transformSpline(*this, spline);
  }

  /*!
   * Applies the operator to a set of coefficients (representing a polynomial on
   * one interval).
   *
   * @param input The polynomial coefficients.
   * @param grid The global grid with respect to which the splines are defined.
   * @param intervalIndex The index of the begin of the interval with respect to
   * the global grid.
   * @tparam T The datatype of the spline.
   * @tparam size The size of the array, i. e. the number of coefficients.
   */
  template <size_t size>
  auto transform(const std::array<T, size> &input, const support::Grid<T> &grid,
                 size_t intervalIndex) const {
    static_assert(size >= 1);
    constexpr size_t OUTPUT_SIZE = outputOrder(size - 1) + 1;

    if (_s.getSupport().getGrid() != grid) {
      throw exceptions::BSplineException(ErrorCode::DIFFERING_GRIDS);
    }

    const auto relativeIndex =
        _s.getSupport().relativeFromAbsolute(intervalIndex);

    auto retVal = internal::make_array<T, OUTPUT_SIZE>(static_cast<T>(0));

    if (relativeIndex) {
      // The interval is part of the Spline's support.
      const auto &coeffs = _s.getCoefficients()[*relativeIndex];

      for (size_t i = 0; i < input.size(); i++) {
        for (size_t j = 0; j < coeffs.size(); j++) {
          retVal[i + j] += input[i] * coeffs[j];
        }
      }
    }
    return retVal;
  }
};
}  // namespace bspline::operators
#endif  // BSPLINE_OPERATORS_SPLINEOPERATOR_H
