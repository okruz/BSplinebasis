#ifndef OKRUZ_BSPLINE_EXAMPLES_MISC_H
#define OKRUZ_BSPLINE_EXAMPLES_MISC_H

#ifdef OKRUZ_EXAMPLES_USE_BOOST_MULTIPRECISION
// Just as a proof of concept.
// Beware boost multiprecision tries to include the eigen headers from
// <Eigen/Core>. Without additional include definitions, they may be placed
// under <eigen3/Eigen/Core>, however, depending on your installation.
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#else
#include <eigen3/Eigen/Dense>
#endif

#include <assert.h>
#include <okruz/bspline/Core.h>

#include <complex>
#include <functional>
#include <vector>

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

namespace okruz::bspline::examples {

#ifdef OKRUZ_EXAMPLES_USE_BOOST_MULTIPRECISION
using data_t = boost::multiprecision::cpp_bin_float_oct;
#else
using data_t = double;
#endif

using DeMat = Eigen::Matrix<data_t, Eigen::Dynamic, Eigen::Dynamic>;

constexpr size_t SPLINE_ORDER = 5;

using Spline = okruz::bspline::Spline<data_t, SPLINE_ORDER>;

struct Eigenspace {
  data_t energy;
  Spline wavefunction;
};

/**
 * @brief setUpSymmetricMatrix Sets up a symmetric matrix, where the matrix
 * elements are defined by a give function.
 * @param b The BilinearForm.
 * @param basis The basis functions (i.e. BSplines).
 * @tparam B The type of the bilinear form.
 * @returns The matrix.
 */
template <typename B>
DeMat setUpSymmetricMatrix(const B &b, const std::vector<Spline> &basis) {
  DeMat ret = DeMat::Zero(basis.size(), basis.size());
  for (size_t i = 0; i < basis.size(); i++) {
    for (size_t j = i; j < basis.size(); j++) {
      const data_t val = b.integrate(basis.at(i), basis.at(j));
      ret(i, j) = val;
      ret(j, i) = val;
    }
  }
  return ret;
}

/**
 * Turns the real Eigen vector into a std::vector.
 */
template <typename Ev>
std::vector<data_t> toStdVector(const Ev &d) {
  std::vector<data_t> ret;
  ret.reserve(d.size());
  for (int i = 0; i < d.size(); i++) {
    ret.push_back(d(i));
  }
  return ret;
}

}  // namespace okruz::bspline::examples
#endif  // OKRUZ_BSPLINE_EXAMPLES_MISC_H
