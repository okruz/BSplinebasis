#ifndef OKRUZ_BSPLINE_EXAMPLES_HYDROGEN_H
#define OKRUZ_BSPLINE_EXAMPLES_HYDROGEN_H

#include <okruz/bspline/Spline.h>

#include <vector>

#include "misc.h"

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

namespace okruz::bspline::examples::hydrogen {

/**
 * The orbital quantum number.
 */
static constexpr int L = 1;

/**
 * Solves the radial hydrogen problem with a BSpline
 * basis.
 * @return Returns the (real) eigenenergies and the corresponding wavefunctions.
 */
std::vector<Eigenspace> solveRadialHydrogen();

}  // namespace okruz::bspline::examples::hydrogen
#endif  // OKRUZ_BSPLINE_EXAMPLES_HYDROGEN_H
