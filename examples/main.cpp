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
#include <fstream>
#include <iostream>
#include <vector>

using namespace okruz::bspline::examples::harmonic_oscillator;

int main() {
  std::cout.precision(20);
  std::vector<HarmonicOscillatorRetVal> harmonicOscillator =
      solveHarmonicOscillator();

  std::cout << "Harmonic Oscillator eigenvalues:\n";
  for (size_t i = 0; i < harmonicOscillator.size(); i++) {
    std::cout << i << "\t" << harmonicOscillator[i].energy << '\n';
  }
  std::cout << '\n' << std::endl;

  {
    std::ofstream out("harmonic-oscillator-wavefunctions.txt");
    out.precision(20);
    for (double x = -10.1; x <= 10.1; x += 0.01) {
      out << x;
      for (const auto &r : harmonicOscillator) {
        out << "\t" << r.wavefunction(x) + r.energy;
      }
      out << "\n";
    }
  }

  return 0;
}
