#ifndef OKRUZ_BSPLINE_OPERATORS_GENERICOPERATOR_H
#define OKRUZ_BSPLINE_OPERATORS_GENERICOPERATOR_H
/*
 * This file contains additional numerical interpolation routines for the
 * splines defined in spline-template.h. The linear algebra routines can
 * be supplied via the implementation of a Solver class (a subclass of
 * internal::ISolver<T>. Implementations based on armadillo and eigen are
 * provided and can be used by defining
 * OKRUZ_BSPLINE_INTERPOLATION_USE_ARMADILLO or
 * OKRUZ_BSPLINE_INTERPOLATION_USE_EIGEN , respectively.
 *
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
using Spline = okruz::bspline::Spline;

/**
 * Marker interface for operators. All proper operators must derive from this
 * interface.
 */
class Operator {};

/**
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
          typename = std::enable_if_t<std::is_base_of<Operator, O>::value>>
decltype(auto) transformSpline(O op, const Spline<T, order> &spline) {
  using OutputArray = decltype(op * spline.getCoefficients().front());

  std::vector<OutputArray> newCoefficients;
  newCoefficients.reserve(spline.getCoefficients().size());

  for (const auto &oldCoeffs : spline.getCoefficients()) {
    newCoefficients.push_back(op * oldCoeffs);
  }

  return Spline(spline.getSupport(), std::move(newCoefficients));
}

}  // namespace okruz::bspline::operators
#endif  // OKRUZ_BSPLINE_OPERATORS_GENERICOPERATOR_H
