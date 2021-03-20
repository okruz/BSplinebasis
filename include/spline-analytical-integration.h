/*
 * This file contains additional integration routines computing a few special integrals for the
 * splines defined in spline-template.h. The integrals are computed from analytical formulas and
 * introduce no additional dependencies.
 *
 * All methods accessing two splines assume that these splines are defined on the same grid (i.e. that
 * both splines have the same interval boundaries within the intersection of their respective supports).
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

#ifndef SPLINE_ANALYTICAL_INTEGRATION_H
#define SPLINE_ANALYTICAL_INTEGRATION_H
#include <spline-template.h>

namespace myspline {
namespace internal {
/*!
 * Efficient integer power for arbitrary type.
 *
 * @param a Basis.
 * @param n Integer exponent.
 * @tparam T Datatype.
 */
template<typename T>
T pow(T a, size_t n) {
    size_t nc = 1;
    T ret = static_cast<T>(1);
    while (nc <= n) {
        if ((nc & n) == nc) {
            ret *= a;
        }
        nc *= 2;
        a *= a;
    }
    return ret;
};


/*!
 * Performs an integral over splines m1 and m2 on one interval. The type of the integral is defined by the integration function f.
 *
 * @param f Function defining the type of integral
 * @param coeffsa First spline's coefficients on the interval of interest.
 * @param coeffsb Second spline's coefficients on the interval of interest.
 * @param x0 Beginning of the interval.
 * @param x1 End of the interval.
 * @tparam T Datatype of both splines.
 * @tparam F Type of Function object f.
 * @tparam sizea Number of coefficients per interval for the first spline.
 * @tparam sizeb Number of coefficients per interval for the second spline.
 */
template<typename T, typename F, size_t sizea, size_t sizeb>
T integrateInterval_analytically(F f, const std::array<T, sizea> &coeffsa, const std::array<T, sizeb> &coeffsb, const T& x0, const T& x1) {
    T result = static_cast<T>(0);;
    const T dxhalf = (x1-x0)/static_cast<T>(2);
    const T xm = (x1 + x0) / static_cast<T>(2);
    for (size_t i = 0; i < sizea; i++) {
        for(size_t j =0;j < sizeb; j++) {
            result += f(i, j,coeffsa[i], coeffsb[j], dxhalf, xm);
        }
    }
    return result;
}; 


/*!
 * Performs an integral over the common support of splines m1 and m2. The type of the integral is defined by the integration function f.
 *
 * @param f Function defining the type of the integral.
 * @param m1 First spline.
 * @param m2 Second spline.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of the first spline.
 * @tparam order2 Order of the second spline.
 */
template<typename T, typename F, size_t order1, size_t order2>
T helper_analytic_integration(F f, const myspline<T, order1> &m1, const myspline<T, order2> &m2){
    size_t startindex1, startindex2, nintervals;
    findOverlappingIntervals(m1, m2, startindex1, startindex2, nintervals);

    if(nintervals == 0) return static_cast<T>(0); // no overlap

    T result = static_cast<T>(0);

    for (size_t interv = 0; interv < nintervals; interv++) {
        result += integrateInterval_analytically<T, F, order1+1, order2 + 1>(f, m1.getCoefficients()[startindex1 + interv],
            m2.getCoefficients()[startindex2+interv],
            m1.getIntervals()[startindex1+interv],
            m1.getIntervals()[startindex1+interv+1]);
    }
    return result;   
};

}; // end namespace internal

/*!
 * Returns the integral \\int\\limits_{-\\infty}^{\\infty} dx m(x). Calculated analytically.
 * 
 * @param m Spline m(x) to be integrated.
 * @tparam T Datatype of the spline m.
 * @tparam order Order of the spline m.
 */
template<typename T, size_t order>
T integrate(const myspline<T, order> &m) {
    T retval = static_cast<T>(0);
    const auto &ints = m.getIntervals();
    for(size_t i = 0; i +1 < ints.size(); i++) {
        const T &start = ints[i];
        const T &end = ints[i+1];
        T pot = (end-start)/static_cast<T>(2); // power of dxhalf, initialised to dxhalf^1
        const T dxhalf_squared = pot * pot;
        const auto &coeffs = m.getCoefficients()[i];
        for(size_t index = 0; index < order + 1; index += 2) {
            retval += static_cast<T>(2) * coeffs[index] * pot / static_cast<T>(index + 1);
            pot *= dxhalf_squared;
        }
    }
    return retval;
};

/*!
 * Returns the integral \\int\\limits_{-\\infty}^{\\infty} dx m1(x) m2(x). Calculated analytically.
 * 
 * @param m1 First spline.
 * @param m2 Second spline.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of the first spline m1.
 * @tparam order2 Order of the second spline m2.
 */
template<typename T, size_t order1, size_t order2>
T overlap(const myspline<T, order1> &m1, const myspline<T, order2> &m2) {
    static constexpr auto f = [](size_t i, size_t j, const T& coeffa, const T& coeffb, const T& dxhalf, [[maybe_unused]] const T& xm) {
        if ((i + j +1) % 2 == 0) return static_cast<T>(0);
        return static_cast<T>(2) * coeffa * coeffb * internal::pow<T>(dxhalf, i + j + 1) / static_cast<T>(i+j+1);
    };
    return internal::helper_analytic_integration(f, m1, m2);
};


/*!
 * Returns the integral \\int\\limits_{-\\infty}^{\\infty} dx m1(x) x m2(x). Calculated analytically.
 * 
 * @param m1 First spline.
 * @param m2 Second spline.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of the first spline m1.
 * @tparam order2 Order of the second spline m2.
 */
template<typename T, size_t order1, size_t order2>
T integrate_x(const myspline<T, order1> &m1, const myspline<T, order2> &m2) {
    static constexpr auto f = [](size_t i, size_t j, const T& coeffa, const T& coeffb, const T& dxhalf, const T& xm) {
        if ((i + j + 1) % 2 == 1) return static_cast<T>(2) * coeffa * coeffb * xm * internal::pow<T>(dxhalf, i +j + 1)/static_cast<T>(i + j +1);
        else return static_cast<T>(2) * coeffa * coeffb * internal::pow<T>(dxhalf, i + j + 2)/static_cast<T>(i+j+2);
    };
    return internal::helper_analytic_integration(f, m1, m2);
};


/*!
 * Returns the integral \\int\\limits_{-\\infty}^{\\infty} dx m1(x) x^2 m2(x). Calculated analytically.
 * 
 * @param m1 First spline.
 * @param m2 Second spline.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of the first spline m1.
 * @tparam order2 Order of the second spline m2.
 */
template<typename T, size_t order1, size_t order2>
T integrate_x2(const myspline<T, order1> &m1, const myspline<T, order2> &m2) {
    static constexpr auto f = [](size_t i, size_t j,const T& coeffa, const T& coeffb, const T& dxhalf, const T& xm) {
        if ((i + j + 2) % 2  == 1) return static_cast<T>(4)*coeffa * coeffb * xm * internal::pow<T>(dxhalf, i + j + 2)/static_cast<T>(i+j+2);
        else return static_cast<T>(2) * coeffa * coeffb * internal::pow<T>(dxhalf,i+j+1) * (dxhalf*dxhalf/static_cast<T>(i+j+3) + xm*xm/ static_cast<T>(i+j+1));
    };
    return internal::helper_analytic_integration(f, m1, m2);
};


/*!
 * Returns the integral \\int\\limits_{-\\infty}^{\\infty} dx m1(x) \frac{\partial}{\partial x} m2(x). Calculated analytically. Assumes m2(x) is continous.
 * 
 * @param m1 First spline.
 * @param m2 Second spline.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of the first spline m1.
 * @tparam order2 Order of the second spline m2.
 */
template<typename T, size_t order1, size_t order2>
T integrate_dx(const myspline<T, order1> &m1, const myspline<T, order2> &m2) {
    static constexpr auto f = [](size_t i, size_t j, const T& coeffa, const T& coeffb, const T& dxhalf, [[maybe_unused]] const T& xm) {
        if (j == 0 || (i+j) % 2 == 0) return static_cast<T>(0);
        else return static_cast<T>(2* j) * coeffa * coeffb * internal::pow<T>(dxhalf, i+j) / static_cast<T>(i+j);
    };
    return internal::helper_analytic_integration(f, m1, m2);
};


/*!
 * Returns the integral \\int\\limits_{-\\infty}^{\\infty} dx x m1(x) \frac{\partial}{\partial x} m2(x). Calculated analytically. Assumes m2(x) is continous.
 * 
 * @param m1 First spline.
 * @param m2 Second spline.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of the first spline m1.
 * @tparam order2 Order of the second spline m2.
 */
template<typename T, size_t order1, size_t order2>
T integrate_x_dx(const myspline<T, order1> &m1, const myspline<T, order2> &m2) {
    static constexpr auto f = [](size_t i, size_t j, const T& coeffa, const T& coeffb, const T& dxhalf, const T& xm) {
        if (j == 0) return static_cast<T>(0);
        else if ((i+j) % 2 == 0) return static_cast<T>(2 * j ) * coeffa * coeffb  * internal::pow<T>(dxhalf, i+j+1) / static_cast<T>(i+j+1);
        else return static_cast<T>(2 * j) * xm * coeffa * coeffb  * internal::pow<T>(dxhalf, i+j) / static_cast<T>(i+j);
    };
    return internal::helper_analytic_integration(f, m1, m2);
};


/*!
 * Returns the integral \\int\\limits_{-\\infty}^{\\infty} dx m1(x) \frac{\partial^2}{\partial x^2} m2(x). Calculated analytically. Assumes m2(x) is at least once continously differentiable.

 * @param m1 First spline.
 * @param m2 Second spline.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of the first spline m1.
 * @tparam order2 Order of the second spline m2.
 */
template<typename T, size_t order1, size_t order2>
T integrate_dx2(const myspline<T, order1> &m1, const myspline<T, order2> &m2) {
    static constexpr auto f = [](size_t i, size_t j, const T& coeffa, const T& coeffb, const T& dxhalf, [[maybe_unused]] const T& xm) {
        if (j < 2 || (i + j) % 2 == 1) return static_cast<T>(0);
        return static_cast<T>(2 * j * (j-1)) * coeffa * coeffb * internal::pow<T>(dxhalf, i+j-1) / static_cast<T>(i+j-1);
    };
    return internal::helper_analytic_integration(f, m1, m2);
}

/*!
 * Returns the integral \\int\\limits_{-\\infty}^{\\infty} dx m1(x) x \frac{\partial^2}{\partial x^2} m2(x). Calculated analytically. Assumes m2(x) is at least once continously differentiable.
 *
 * @param m1 First spline.
 * @param m2 Second spline.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of the first spline m1.
 * @tparam order2 Order of the second spline m2.
 */
template<typename T, size_t order1, size_t order2>
T integrate_x_dx2(const myspline<T, order1> &m1, const myspline<T, order2> &m2) {
    static constexpr auto f = [](size_t i, size_t j, const T& coeffa, const T& coeffb, const T& dxhalf, const T& xm) {
        if (j < 2) return static_cast<T>(0);
        else if ((i+j) % 2 == 1) return static_cast<T>(2 * j * (j-1)) * coeffa * coeffb * internal::pow<T>(dxhalf, i + j) / static_cast<T>(i+j);
        else return static_cast<T>(2 * j * (j-1)) * coeffa * coeffb * xm *  internal::pow<T>(dxhalf, i + j -1)/static_cast<T>(i + j -1);
    };
    return internal::helper_analytic_integration(f, m1, m2);
};


/*!
 * Returns the integral \\int\\limits_{-\\infty}^{\\infty} dx m1(x) x^2 \frac{\partial^2}{\partial x^2} m2(x). Calculated analytically. Assumes m2(x) is at least once continously differentiable.
 *
 * @param m1 First spline.
 * @param m2 Second spline.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of the first spline m1.
 * @tparam order2 Order of the second spline m2.
 */
template<typename T, size_t order1, size_t order2>
T integrate_x2_dx2(const myspline<T, order1> &m1, const myspline<T, order2> &m2) {
    static constexpr auto f = [](size_t i, size_t j, const T& coeffa, const T& coeffb, const T& dxhalf, const T& xm) {
        if (j < 2) return static_cast<T>(0);
        else if ((i+j) % 2 == 1) return static_cast<T>(4 * j * (j-1)) * xm * coeffa * coeffb * internal::pow<T>(dxhalf, i + j) / static_cast<T>(i+j);
        else return static_cast<T>(2 * j * (j-1)) * coeffa * coeffb * internal::pow<T>(dxhalf, i + j - 1) * (dxhalf * dxhalf /static_cast<T>(i + j +1) + xm * xm /static_cast<T>(i + j -1));
    };
    return internal::helper_analytic_integration(f, m1, m2);
};

}; // end  namespace myspline
#endif // SPLINE_ANALYTICAL_INTEGRATION_H
