#ifndef OKRUZ_BSPLINE_INTEGRATION_NUMERICAL_H
#define OKRUZ_BSPLINE_INTEGRATION_NUMERICAL_H
/*
 * This file contains an additional numerical integration routine for the
 * Splines based on the Gauss-Legendre routines in boost/math/quadrature.
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

#include <boost/math/quadrature/gauss.hpp>
#include <okruz/bspline/Spline.h>
#include <okruz/bspline/internal/misc.h>

namespace okruz::bspline::integration {
using namespace boost::math::quadrature;

/*!
 * Calculates the 1D integral \\int\\limits_{-\\infty}^{\\infty} dx m1(x) f(x)
 * m2(x). Calcualted numerically with each interval of the common support of the
 * two splines being integrated separately using a Gauss-Legendre scheme of
 * order ordergl. Uses the boost Gauss-Legendre fixed point integration scheme.
 *
 * @param f Function f(x) to be multiplied to the integrand.
 * @param m1 First spline m1(x).
 * @param m2 Second spline m2(x).
 * @tparam ordergl Order of the Gauss-Legendre integration scheme provided by
 * the boost library.
 * @tparam T Datatype of the calculation.
 * @tparam order1 Order of the spline m1(x).
 * @tparam order2 Order of the spline m2(x).
 */
template <size_t ordergl, typename T, typename F, size_t order1, size_t order2>
T integrate(const F &f, const okruz::bspline::Spline<T, order1> &m1,
            const okruz::bspline::Spline<T, order2> &m2) {
  using okruz::bspline::support::Support;
  Support newSupport = m1.getSupport().calcIntersection(m2.getSupport());
  const size_t nintervals = newSupport.numberOfIntervals();
  if (nintervals == 0)
    return static_cast<T>(0); // no overlap

  T result = static_cast<T>(0);
  for (size_t interv = 0; interv < nintervals; interv++) {

    const auto ai = newSupport.absoluteFromRelative(interv);
    const auto m1Index = m1.getSupport().relativeFromAbsolute(ai).value();
    const auto m2Index = m2.getSupport().relativeFromAbsolute(ai).value();

    const T &xstart = m1.getSupport().at(m1Index);
    const T &xend = m1.getSupport().at(m1Index + 1);
    const T xm = (xstart + xend) / static_cast<T>(2);
    const auto &c1 = m1.getCoefficients().at(m1Index);
    const auto &c2 = m2.getCoefficients().at(m2Index);
    const auto fwrap = [&c1, &c2, &xm, &f](const T &x) {
      using namespace okruz::bspline::internal;
      return f(x) * evaluateInterval(x, c1, xm) * evaluateInterval(x, c2, xm);
    };
    result += gauss<T, ordergl>::integrate(fwrap, xstart, xend);
  }
  return result;
};
};     // namespace okruz::bspline::integration
#endif // OKRUZ_BSPLINE_INTEGRATION_NUMERICAL_H
