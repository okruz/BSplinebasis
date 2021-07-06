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

#include "harmonic-oscillator.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <okruz/bspline/BSplineGenerator.h>
#include <okruz/bspline/integration/analytical.h>

#include <algorithm>
#include <assert.h>
#include <complex>
#include <functional>

namespace okruz::bspline::examples::harmonic_oscillator {

using namespace okruz::bspline;
using DeMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * @brief setUpKnotsVector Sets up the knots vector on which the basis splines
 * are defined.
 * @return A vector of doubles representing the knots.
 */
static std::vector<double> setUpKnotsVector() {
  std::vector<double> ret{0};

  for (int i = 1; i < 100; i++) {
    const double val = double(i * i) / 1000.0;
    ret.push_back(val);
    ret.push_back(-val);
  }
  std::sort(ret.begin(), ret.end());
  return ret;
}

/**
 * @brief setUpBasis Sets up the BSpline basis.
 * @return A vector of BSplines representing the basis.
 */
static std::vector<Spline> setUpBasis() {
  BSplineGenerator gen(setUpKnotsVector());
  return gen.template generateBSplines<SPLINE_ORDER + 1>();
}

/**
 * @brief setUpSymmetricMatrix Sets up a symmetric matrix, where the matrix
 * elements are defined by a give function.
 * @param f The function defining the matrix elements.
 * @param basis The basis functions (i.e. BSplines).
 * @return The matrix.
 */
static DeMat setUpSymmetricMatrix(
    const std::function<double(const Spline &, const Spline &)> &f,
    const std::vector<Spline> &basis) {
  DeMat ret = DeMat::Zero(basis.size(), basis.size());
  for (size_t i = 0; i < basis.size(); i++) {
    for (size_t j = i; j < basis.size(); j++) {
      double val = f(basis.at(i), basis.at(j));
      ret(i, j) = val;
      ret(j, i) = val;
    }
  }
  return ret;
}

/**
 * @brief getIdentityPerm Returns the identity permutation for size elements.
 * @param size The number of elements.
 * @return Returns the identity permutation for size elements.
 */
static std::vector<size_t> getIdentityPerm(size_t size) {
  std::vector<size_t> ret(size);
  for (size_t i = 0; i < size; i++)
    ret[i] = i;
  return ret;
}

/**
 * Turns the complex Eigen vector into a real std::vector.
 */
template <typename Ev> static std::vector<double> toRealVector(const Ev &d) {
  std::vector<double> ret;
  ret.reserve(d.size());
  for (int i = 0; i < d.size(); i++) {
    ret.push_back(d(i).real());

    // Imaginary part should be negligible.
    assert(abs(d(i).imag()) < 1.0e-10);
  }
  return ret;
}

std::vector<HarmonicOscillatorRetVal> solveHarmonicOscillator() {
  // Get the basis.
  std::vector<Spline> basis = setUpBasis();

  // kinetic term -1/2 d^2/dx^2
  DeMat hamiltonian =
      -0.5 * setUpSymmetricMatrix(okruz::bspline::integration::integrate_dx2<
                                      double, SPLINE_ORDER, SPLINE_ORDER>,
                                  basis);

  // potential term 1/2 x^2
  hamiltonian +=
      0.5 * setUpSymmetricMatrix(
                okruz::bspline::integration::integrate_x2<double, SPLINE_ORDER,
                                                          SPLINE_ORDER>,
                basis);

  // overlap matrix
  DeMat overlapMatrix = setUpSymmetricMatrix(
      okruz::bspline::integration::overlap<double, SPLINE_ORDER, SPLINE_ORDER>,
      basis);

  // Solve the generalized eigenvalue problem A.x = lambda B.x
  Eigen::GeneralizedEigenSolver<DeMat> ges;
  ges.compute(hamiltonian, overlapMatrix);

  // Retrieve the (complex) eigenvalues and eigenvectors.
  // As the overlap matrix is symmetric and positively definite, the imaginary
  // parts of both should be zero within numerical accuracy.
  const auto eigenvalues = ges.eigenvalues();
  const auto eigenvectors = ges.eigenvectors();

  // Get the identity permutation.
  std::vector<size_t> perm = getIdentityPerm(basis.size());

  // Get the sorted permutation (by the real part of the eigenvalues).s
  std::sort(perm.begin(), perm.end(), [&eigenvalues](size_t i, size_t j) {
    const double vali = eigenvalues(i).real();
    const double valj = eigenvalues(j).real();
    return vali < valj;
  });

  std::vector<HarmonicOscillatorRetVal> ret;
  ret.reserve(10);
  // Return the eienvalues and eigenfunctions correspondig to the ten lowest
  // eigenvalues.  The imaginary part should be zero within numerical accuracy.
  for (size_t i = 0; i < 10; i++) {
    size_t index = perm[i];
    const auto eigenvalue = eigenvalues(index);
    const auto eigenvector = toRealVector(eigenvectors.col(index));
    auto wavefunction = okruz::bspline::linearCombination(
        eigenvector.begin(), eigenvector.end(), basis.begin(), basis.end());

    // Imaginary part should be negligible.
    assert(abs(eigenvalue.imag()) < 1.0e-10);
    ret.push_back({eigenvalue.real(), std::move(wavefunction)});
  }
  return ret;
}

} // namespace okruz::bspline::examples::harmonic_oscillator
