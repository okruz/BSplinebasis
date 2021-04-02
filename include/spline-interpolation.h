/*
 * This file contains additional numerical interpolation routines for the
 * splines defined in spline-template.h. The numerical routines are based on the
 * linear algebra library Eigen and are relegated to this files to
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

#ifndef SPLINE_INTERPOLATION_H
#define SPLINE_INTERPOLATION_H
#include <spline-template.h>
#include <Eigen/Dense>
#include <array>
namespace myspline {

/*!
 * Represents either the first or last node of the interpolation grid. Used to define boundary conditions.
 */
enum class Node {FIRST, LAST};

/*!
 * Represents a boundary condition, i.e. one fixed derivative on either the first or last node of the interpolation grid.
 */
template <typename T>
struct boundary { 
    Node node = Node::FIRST;        /*! Node to apply the boundary condition to. */ 
    size_t derivative = 1;          /*! Order of the derivative to fix. */
    T value = static_cast<T>(0);    /*! Value of the derivative. */
};

namespace internal {
/*!
 * Generates the default boundary conditions by setting as many derivatives to zero as needed, starting from the first derivative.
 *
 * @tparam T Datatype of the spline.
 * @tparam order Polynomial order of the spline.
 */
template<typename T, size_t order>
std::array<boundary<T>, order-1> defaultBoundaries() {
     static_assert(order >= 1, "Order may not be zero.");
     std::array<boundary<T>,order-1> ret;
     for (size_t i = 0; i < order -1; i++) {
          if (i % 2 == 0) ret[i] = boundary<T>{.node = Node::FIRST, .derivative = i/2 + 1, .value = static_cast<T>(0)};
          else ret[i] = boundary<T>{.node = Node::LAST, .derivative = (i-1)/2 + 1, .value = static_cast<T>(0)};
     }
     return ret;
}

/*!
 * Computes exponent! / (exponent - deriv)!.
 *
 * @param exponent Exponent of the monomial.
 * @param deriv Order of the derivative.
 */
template<typename T>
T faculty_ratio(size_t exponent, size_t deriv) {
    assert(deriv <= exponent);
    T ret = static_cast<T>(exponent);
    for (size_t j = 1; j < deriv; j++) ret *= static_cast<T>(exponent -j);
    return ret;
}
}; // end namespace internal

/*!
 * Interpolates the data given by x and y with a spline of order order. order-1 additional conditions are needed for a well defined problem. These can be supplied by fixing derivatives on the first and last node.
 * 
 * @param x Data on the abscissa. The grid points must be in (steadily) increasing order.
 * @param y Data on the ordinate.
 * @param boundaries Boundary conditions.
 * @tparam T Datatype of the spline and data.
 * @tparam order Order of the spline.
 */
template<typename T, size_t order>
myspline<T, order> interpolate(const std::vector<T> &x, const std::vector<T> &y, 
                             const std::array<boundary<T>, order-1> boundaries = internal::defaultBoundaries<T,order>()) {
    static_assert(order >= 1, "Order may not be zero.");
    assert(x.size() >= 2 && x.size() == y.size());
    assert(internal::isSteadilyIncreasing(x));
    using DeMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using DeVec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    constexpr size_t NUM_COEFFS = order +1;    

    DeMat m = DeMat::Zero(NUM_COEFFS*(x.size()-1), NUM_COEFFS*(x.size()-1));
    DeVec b = DeVec::Zero(NUM_COEFFS*(x.size()-1));
    size_t row_counter = 0;
    {
        const T dx1 = (x[0] - x[1])/static_cast<T>(2);
        {
            T power_of_dx1 = static_cast<T>(1);
            for (size_t i = 0; i <= order; i++) {
                m(row_counter, i) = power_of_dx1;
                power_of_dx1 *= dx1;
            }
            b(row_counter) = y.front();
            row_counter++;
        }

        for (const auto &bo: boundaries) {
            assert(bo.derivative <= order);
            if (bo.node == Node::FIRST) {
                T power_of_dx1 = static_cast<T>(1); 
                for (size_t i = bo.derivative; i <= order; i++) {
                    m(row_counter, i) = internal::faculty_ratio<T>(i, bo.derivative) * power_of_dx1;
                    power_of_dx1 *= dx1;
                }
                b(row_counter) =bo.value;
                row_counter++;
            }
        }
    }

    for (size_t c = 1; c + 1 < x.size(); c++) {
        const T dx1 = (x[c] - x[c-1])/static_cast<T>(2);
        const T dx2 = (x[c] - x[c+1])/static_cast<T>(2);

        {
            T power_of_dx1 = static_cast<T>(1); 
            for (size_t i = 0; i <= order; i++) {
                  m(row_counter, NUM_COEFFS * (c-1) + i) = power_of_dx1;
                  power_of_dx1 *= dx1;
            }
            b(row_counter) = y[c];
            row_counter++;
        }

        {
            T power_of_dx2 = static_cast<T>(1); 
            for (size_t i = 0; i <= order; i++) {
                m(row_counter, NUM_COEFFS * c + i) = power_of_dx2;
                power_of_dx2 *= dx2;
            }
            b(row_counter) = y[c];
            row_counter++;
        }

        for(size_t deriv = 1; deriv < order; deriv++) {
            T power_of_dx1 = static_cast<T>(1);
            T power_of_dx2 = static_cast<T>(1);
            for (size_t i = deriv; i <= order; i++) {
                m(row_counter, NUM_COEFFS * (c-1) + i) = internal::faculty_ratio<T>(i, deriv) * power_of_dx1;
                m(row_counter, NUM_COEFFS * c + i) = -internal::faculty_ratio<T>(i, deriv) * power_of_dx2;
                power_of_dx1 *= dx1;
                power_of_dx2 *= dx2;
            }
            row_counter++;
        }
    }

    {
        const T dx2 =  (x.back() - x[x.size()-2])/static_cast<T>(2);
        {
            T power_of_dx2 = static_cast<T>(1);
            for (size_t i = 0; i <= order; i++) {
                m(row_counter, NUM_COEFFS * (x.size()-2) + i) = power_of_dx2;
                power_of_dx2 *= dx2;
            }
            b(row_counter) = y.back();
            row_counter++;
        }

        for (const auto &bo: boundaries) {
            if (bo.node == Node::LAST) {
                T power_of_dx2 = static_cast<T>(1);
                for (size_t i = bo.derivative; i <= order; i++) {
                    m(row_counter, NUM_COEFFS * (x.size() -2) + i) = internal::faculty_ratio<T>(i, bo.derivative) * power_of_dx2;
                    power_of_dx2 *= dx2;
                }
                b(row_counter) = bo.value;
                row_counter++;
            }
        }

    }

    assert(row_counter == NUM_COEFFS * (x.size() -1));

    DeVec result = m.colPivHouseholderQr().solve(b);
    std::vector<std::array<T, NUM_COEFFS>> coeffs((x.size() -1));
    for (size_t i = 0 ; i +1 < x.size(); i++) {
        std::array<T, NUM_COEFFS> &coeffsi = coeffs[i];
        for(size_t j = 0; j < NUM_COEFFS; j++) coeffsi[j] = result(NUM_COEFFS * i + j);
    }
    return myspline<T, order>(x, std::move(coeffs));
}
};
#endif //SPLINE_INTERPOLATION_H
