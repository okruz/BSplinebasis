/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include "spline-potential.h"

#include <bspline/BSplineGenerator.h>
#define BSPLINE_INTERPOLATION_USE_EIGEN
#include <bspline/interpolation/interpolation.h>

#include <algorithm>
#include <eigen3/Eigen/Eigenvalues>

namespace bspline::examples::spline_potential {

using namespace bspline;

Spline interpolateFunction(std::vector<data_t> gridPoints,
                           const std::function<data_t(data_t)> &func) {
  support::Support<data_t> support{Grid<data_t>{std::move(gridPoints)},
                                   support::Construction::WHOLE_GRID};
  std::vector<data_t> y;
  y.reserve(support.size());
  for (const auto x : support) {
    y.push_back(func(x));
  }
  return interpolation::interpolate_using_eigen<data_t, SPLINE_ORDER>(
      std::move(support), y);
}

/**
 * @brief setUpBasis Sets up the BSpline basis.
 * @param[in] grid The global grid on which to generate the splines.
 * @return A vector of BSplines representing the basis.
 */
static std::vector<Spline> setUpBasis(const support::Grid<data_t> &grid) {
  const std::vector<data_t> knotsVector{grid.begin(), grid.end()};
  BSplineGenerator gen(knotsVector, grid);
  return gen.template generateBSplines<SPLINE_ORDER + 1>();
}

std::vector<Eigenspace> solveSEWithSplinePotential(Spline v) {
  // Get the basis.
  const std::vector<Spline> basis = setUpBasis(v.getSupport().getGrid());

  // Hamiltonian operator -1/2 d^2/dx^2 + v(x)
  auto hamiltonOperator = (static_cast<data_t>(-1) / 2) * operators::Dx<2>{} +
                          operators::SplineOperator{std::move(v)};

  const DeMat hamiltonian = setUpSymmetricMatrix(
      integration::BilinearForm{std::move(hamiltonOperator)}, basis);

  // overlap matrix
  const DeMat overlapMatrix =
      setUpSymmetricMatrix(integration::ScalarProduct{}, basis);

  // Solve the generalized eigenvalue problem A.x = lambda B.x
  const Eigen::GeneralizedSelfAdjointEigenSolver<DeMat> ges{hamiltonian,
                                                            overlapMatrix};

  // Retrieve the (real) eigenvalues and eigenvectors.
  const auto &eigenvalues = ges.eigenvalues();
  const auto &eigenvectors = ges.eigenvectors();

  std::vector<Eigenspace> ret;
  ret.reserve(10);
  // Return the eienvalues and eigenfunctions correspondig to the ten lowest
  // eigenvalues.
  for (size_t i = 0; i < 10; i++) {
    const auto eigenvalue = eigenvalues(i);
    const auto eigenvector = toStdVector(eigenvectors.col(i));
    auto wavefunction = bspline::linearCombination(eigenvector, basis);

    ret.push_back({eigenvalue, std::move(wavefunction)});
  }
  return ret;
}

}  // namespace bspline::examples::spline_potential
