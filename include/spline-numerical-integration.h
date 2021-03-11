/*
 * This file contains additional numerical integration routines (1D, 2D and 3D) for the
 * splines defined in spline-template.h. The numerical routines are based on the
 * Gauss-Legendre routines in boost/math/quadrature and are relegated to this files to
 * minimise the dependencies of spline-template.h.
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

#ifndef SPLINE_NUMERICAL_INTEGRATION_H
#define SPLINE_NUMERICAL_INTEGRATION_H
#include <spline-template.h>
#include <boost/math/quadrature/gauss.hpp> 

namespace myspline {
using namespace boost::math::quadrature;

namespace internal {
/*
 * Normally, for every evaluation of a spline, a binary search for the correct interval is necessary.
 * This method is defined in order to integrate every interval separately during the 1D integration, omitting the necessity for the binary search.
 */
template<typename T, size_t ARRAY_SIZE>
T evaluateInterval(const T& x, const std::array<T, ARRAY_SIZE> &coeffs, const T& xm) {
    T result = static_cast<T>(0), xpot = static_cast<T>(1);
    const T dx = x - xm;
    for (const T& c: coeffs) {
        result += c * xpot;
        xpot *= dx;
    }
    return result;
};
}; //end namespace internal


/*!
 * Calculates the 1D integral \\int\\limits_{-\\infty}^{\\infty} dx m1(x) f(x) m2(x).
 * Calcualted numerically with each interval of the common support of the two splines being integrated separately using a Gauss-Legendre scheme of order ordergl.
 * Uses the boost Gauss-Legendre fixed point integration scheme.
 * 
 * @param f Function f(x) to be multiplied to the integrand.
 * @param m1 First spline m1(x).
 * @param m2 Second spline m2(x).
 * @tparam T Datatype of the calculation.
 * @tparam order1 Order of the spline m1(x).
 * @tparam order2 Order of the spline m2(x).
 * @tparam ordergl Order of the Gauss-Legendre integration scheme provided by the boost library.
 */
template<typename T, size_t order1, size_t order2, size_t ordergl>
T integrate(const std::function<T(const T&)> &f, const myspline<T, order1> &m1, const myspline<T, order2> &m2) {
    size_t startindex1, startindex2, nintervals;
    internal::findOverlappingIntervals(m1, m2, startindex1, startindex2, nintervals);
    if (nintervals == 0) return static_cast<T>(0); // no overlap
        
    T result = static_cast<T>(0);
    for (size_t interv = 0; interv < nintervals; interv++) { // Integrate every interval on its own as Gauss-Legendre-integration is exact for polynomials of degree 2*neval -1 or less
        const T& xstart = m1.getIntervals().at(startindex1 + interv);
        const T& xend = m1.getIntervals().at(startindex1 + interv +1);
        const T xm = (xstart + xend) / static_cast<T>(2);
        const auto &c1 = m1.getCoefficients().at(startindex1 + interv);
        const auto &c2 = m2.getCoefficients().at(startindex2 + interv);
        const auto fwrap = [&c1, &c2, &xm, &f](const T &x) {
            return f(x) * internal::evaluateInterval(x, c1, xm) * internal::evaluateInterval(x, c2, xm);
        };
        result += gauss<T, ordergl>::integrate(fwrap, xstart, xend);
    }
    return result;
};
};
#endif //SPLINE_NUMERICAL_INTEGRATION_H
