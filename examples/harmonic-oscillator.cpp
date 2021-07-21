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

#include <okruz/bspline/BSplineGenerator.h>
#include <okruz/bspline/integration/analytical.h>

#include <algorithm>
#include <eigen3/Eigen/Eigenvalues>

namespace okruz::bspline::examples::harmonic_oscillator {

using namespace okruz::bspline;

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

std::vector<Eigenspace> solveHarmonicOscillator() {
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
  Eigen::GeneralizedSelfAdjointEigenSolver<DeMat> ges;
  ges.compute(hamiltonian, overlapMatrix);

  // Retrieve the (real) eigenvalues and eigenvectors.
  const auto eigenvalues = ges.eigenvalues();
  const auto eigenvectors = ges.eigenvectors();

  // Get the identity permutation.
  std::vector<size_t> perm = getIdentityPerm(basis.size());

  // Get the sorted permutation (by the value of the corresponding eigenvalues).
  std::sort(perm.begin(), perm.end(), [&eigenvalues](size_t i, size_t j) {
    return eigenvalues(i) < eigenvalues(j);
  });

  std::vector<Eigenspace> ret;
  ret.reserve(10);
  // Return the eienvalues and eigenfunctions correspondig to the ten lowest
  // eigenvalues.
  for (size_t i = 0; i < 10; i++) {
    size_t index = perm[i];
    const auto eigenvalue = eigenvalues(index);
    const auto eigenvector = toStdVector(eigenvectors.col(index));
    auto wavefunction = okruz::bspline::linearCombination(
        eigenvector.begin(), eigenvector.end(), basis.begin(), basis.end());

    ret.push_back({eigenvalue, std::move(wavefunction)});
  }
  return ret;
}

}  // namespace okruz::bspline::examples::harmonic_oscillator
