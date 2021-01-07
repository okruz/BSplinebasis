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
            return f(x) * evaluateInterval(x, c1, xm) * evaluateInterval(x, c2, xm);
        };
        result += gauss<T, ordergl>::integrate(fwrap, xstart, xend);
    }
    return result;
};


/*!
 * Calculates the integral \\int\\limits_{-\\infty}^{\\infty} dx \\int\\limits_{-\\infty}^{\\infty} dy m1x(x) m1y(y) f(x, y) m2x(x) m2y(y).
 * Calculated numerically with ~(ordergl-1)/2 evaluations per dimension, distributed over the support of the splines.
 * Uses the boost Gauss-Legendre fixed point integration scheme.
 * 
 * @param f Function f(x, y) to be multiplied to the integrand.
 * @param m1x First spline m1(x) for the x coordinate.
 * @param m2x Second spline m2(x) for the x coordinate.
 * @param m1y First spline m1(y) for the y coordinate.
 * @param m2y Second spline m2(y) for the y coordinate.
 * @tparam T Datatype of the calculation.
 * @tparam order1x Order of the spline m1(x).
 * @tparam order2x Order of the spline m2(x).
 * @tparam order1y Order of the spline m1(y).
 * @tparam order2y Order of the spline m2(y).
 * @tparam ordergl Order of the (1D) Gauss-Legendre integration scheme provided by the boost library.
 * @deprecated This method may be removed in the future.
 */
[[deprecated]]
template<typename T, size_t order1x, size_t order2x, size_t order1y, size_t order2y, size_t ordergl>
T integrate2d(const std::function<T(const T&, const T&)> &f, const myspline<T, order1x> &m1x, const myspline<T, order2x> &m2x, const myspline<T, order1y> &m1y, const myspline<T, order2y> &m2y) {
    const T ax = std::max(m1x.start(), m2x.start());
    const T bx = std::min(m1x.end(), m2x.end());
    const T ay = std::max(m1y.start(), m2y.start());
    const T by = std::min(m1y.end(), m2y.end());
    if (ax >= bx || ay >= by) return static_cast<T>(0); // No overlap
    const auto xintegration = [&m1x, &m2x, &m1y, &m2y, &f, &ay, &by] (const T&x) {
        const auto yintegration = [&m1y, &m2y, &f, &x](const T& y) {
            return m1y(y) * m2y(y) *f(x, y);
        };
        return m1x(x) * m2x(x) * gauss<T, ordergl>::integrate(yintegration, ay, by);
    };

    return gauss<T, ordergl>::integrate(xintegration, ax, bx);
};

/*!
 * Calculates the integral \\int\\limits_{-\\infty}^{\\infty} dx \\int\\limits_{-\\infty}^{\\infty} dy \\int\\limits_{-\\infty}^{\\infty} dz m1x(x) m1y(y) m1z(z) f(x, y, z) m2x(x) m2y(y) m2z(z).
 * Calculated numerically with ~(ordergl-1)/2 evaluations per dimension, distributed over the support of the splines.
 * Uses the boost Gauss-Legendre fixed point integration scheme.
 * 
 * @param f Function f(x, y, z) to be multiplied to the integrand.
 * @param m1x First spline m1(x) for the x coordinate.
 * @param m2x Second spline m2(x) for the x coordinate.
 * @param m1y First spline m1(y) for the y coordinate.
 * @param m2y Second spline m2(y) for the y coordinate.
 * @param m1z First spline m1(z) for the z coordinate.
 * @param m2z Second spline m2(x) for the z coordinate.
 * @tparam T Datatype of the calculation.
 * @tparam order1x Order of the spline m1(x).
 * @tparam order2x Order of the spline m2(x).
 * @tparam order1y Order of the spline m1(y).
 * @tparam order2y Order of the spline m2(y).
 * @tparam order1z Order of the spline m1(z).
 * @tparam order2z Order of the spline m2(z).
 * @tparam ordergl Order of the (1D) Gauss-Legendre integration scheme provided by the boost library.
 * @deprecated This method may be removed in the future.
 */
[[deprecated]]
template<typename T, size_t order1x, size_t order2x, size_t order1y, size_t order2y, size_t order1z, size_t order2z, size_t ordergl>
T integrate3d(const std::function<T(const T&, const T&, const T&)> &f, const myspline<T, order1x> &m1x, const myspline<T, order2x> &m2x, const myspline<T, order1y> &m1y, const myspline<T, order2y> &m2y , const myspline<T, order1z> &m1z, const myspline<T, order2z> &m2z) {
    T ax = std::max(m1x.start(), m2x.start());
    T bx = std::min(m1x.end(), m2x.end());
    T ay = std::max(m1y.start(), m2y.start());
    T by = std::min(m1y.end(), m2y.end());
    T az = std::max(m1z.start(), m2z.start());
    T bz = std::min(m1z.end(), m2z.end());
    if (ax >= bx || ay >= by || az >= bz) return static_cast<T>(0); // No overlap
    const auto xintegration = [&m1x, &m2x, &m1y, &m2y, &m1z, &m2z, &f, &ay, &by, &az, &bz] (const T&x) {
        const auto yintegration = [&m1y, &m2y, &m1z, &m2z, &f, &x, &az, &bz](const T& y) {
            const auto  zintegration = [&m1z, &m2z, &f, &x, &y] (const T& z) {
                return m1z(z) * m2z(z) * f(x, y, z);
            };
            return m1y(y) * m2y(y) * gauss<T,ordergl>::integrate(zintegration, az, bz);
        };
        return m1x(x) * m2x(x) * gauss<T, ordergl>::integrate(yintegration, ay, by);
    };
    return gauss<T, ordergl>::integrate(xintegration, ax, bx);
};
};
#endif //SPLINE_NUMERICAL_INTEGRATION_H
