#ifndef OKRUZ_BSPLINE_EXAMPLES_HARMONICOSCILLATOR_H
#define OKRUZ_BSPLINE_EXAMPLES_HARMONICOSCILLATOR_H

#include <okruz/bspline/Spline.h>

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

namespace okruz::bspline::examples::harmonic_oscillator {

constexpr size_t SPLINE_ORDER = 5;

using Spline = okruz::bspline::Spline<double, SPLINE_ORDER>;

struct HarmonicOscillatorRetVal {
  double energy;
  Spline wavefunction;
};

/**
 * Solves the quantum mechanical harmonic oscillator problem with a BSpline
 * basis.
 * @return Returns the (real) eigenenergies and the corresponding wavefunctions.
 */
std::vector<HarmonicOscillatorRetVal> solveHarmonicOscillator();

}  // namespace okruz::bspline::examples::harmonic_oscillator
#endif  // OKRUZ_BSPLINE_EXAMPLES_HARMONICOSCILLATOR_H
