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

#include "hydrogen.h"

#include <okruz/bspline/BSplineGenerator.h>
#include <okruz/bspline/integration/analytical.h>

#include <algorithm>
#include <eigen3/Eigen/Eigenvalues>

namespace okruz::bspline::examples::hydrogen {

using namespace okruz::bspline;

/**
 * @brief setUpKnotsVector Sets up the knots vector on which the basis splines
 * are defined.
 * @return A vector of doubles representing the knots.
 */
static std::vector<double> setUpKnotsVector() {
  std::vector<double> ret;

  size_t numberOfZeros = 1;

  if (SPLINE_ORDER + 1 > L) {
    numberOfZeros = SPLINE_ORDER + 1 - L;
  }

  // Adding a knot multiple times alters the continuity properties of the
  // generated splines at the corresponding grid point (see literature on
  // BSplines). This guarantees that the radial wavefunctions have the correct
  // scaling in the vicinity of r=0 (at least for sufficiently low L).
  for (size_t i = 0; i < numberOfZeros; i++) ret.push_back(0.0);

  // First point of the logarithmic grid.
  const double rmin = 0.01;

  // Last point of the logarithmic grid.
  const double rmax = 1.0e3;

  // Roughly the number of grid points on the logarithmic grid.
  const int numberOfGridPoints = 200;

  // logarithmic step
  const double step = pow(rmax / rmin, 1.0 / double(numberOfGridPoints));

  for (int i = 1; i <= numberOfGridPoints; i++) {
    ret.push_back(rmin * pow(step, i));
  }
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

std::vector<Eigenspace> solveRadialHydrogen() {
  static_assert(L >= 0, "L may not be below zero.");

  // Get the basis.
  std::vector<Spline> basis = setUpBasis();

  // First part of the kinetic term -d^2/dr^2 . Includes the term r^2 from the
  // functional determinant.
  DeMat hamiltonian = -setUpSymmetricMatrix(
      okruz::bspline::integration::integrate_x2_dx2<double, SPLINE_ORDER,
                                                    SPLINE_ORDER>,
      basis);

  // Second part of the kinetic term -2/r d/dr . Includes the term r^2 from the
  // functional determinant.
  hamiltonian +=
      -2 * setUpSymmetricMatrix(
               okruz::bspline::integration::integrate_x_dx<double, SPLINE_ORDER,
                                                           SPLINE_ORDER>,
               basis);

  if constexpr (L != 0) {
    // Third part of the kinetic term L * (L + 1) / r^2 . Includes the term r^2
    // from the functional determinant.
    hamiltonian +=
        L * (L + 1) *
        setUpSymmetricMatrix(
            okruz::bspline::integration::overlap<double, SPLINE_ORDER,
                                                 SPLINE_ORDER>,
            basis);
  }

  // potential term -2/r. Includes the term r^2 from the functional determinant.
  hamiltonian +=
      -2 * setUpSymmetricMatrix(
               okruz::bspline::integration::integrate_x<double, SPLINE_ORDER,
                                                        SPLINE_ORDER>,
               basis);

  // Overlap matrix. Includes the term r^2 from the functional determinant.
  DeMat overlapMatrix = setUpSymmetricMatrix(
      okruz::bspline::integration::integrate_x2<double, SPLINE_ORDER,
                                                SPLINE_ORDER>,
      basis);

  // Solve the generalized eigenvalue problem A.x = lambda B.x
  Eigen::GeneralizedSelfAdjointEigenSolver<DeMat> ges;
  ges.compute(hamiltonian, overlapMatrix);

  // Retrieve the eigenvalues and eigenvectors.
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
    auto wavefunction = okruz::bspline::linearCombination(eigenvector, basis);
    ret.push_back({eigenvalue, std::move(wavefunction)});
  }
  return ret;
}

}  // namespace okruz::bspline::examples::hydrogen
