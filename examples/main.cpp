/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "harmonic-oscillator-spline.h"
#include "harmonic-oscillator.h"
#include "hydrogen.h"

using namespace bspline::examples::hydrogen;
using namespace bspline::examples::harmonic_oscillator;
using namespace bspline::examples;

void harmonicOscillator(bool useSpline) {
  std::cout.precision(std::numeric_limits<data_t>::max_digits10);
  std::vector<Eigenspace> harmonicOscillator;

  if (useSpline) {
    harmonicOscillator = harmonic_oscillator_spline::solveHarmonicOscillator();

  } else {
    harmonicOscillator = solveHarmonicOscillator();
  }

  std::cout << "Harmonic Oscillator eigenvalues:\n\n";
  std::cout << std::setw(2) << "n"
            << "\t" << std::setw(20) << "energy"
            << "\t" << std::setw(20) << "relative deviation" << '\n';
  for (size_t i = 0; i < harmonicOscillator.size(); i++) {
    const data_t expected = static_cast<data_t>(2 * i + 1) / 2;
    const data_t relativeDev =
        abs((harmonicOscillator[i].energy - expected) / expected);
    std::cout << std::setw(2) << i << "\t" << std::setw(20)
              << harmonicOscillator[i].energy << "\t" << std::setw(20)
              << relativeDev << '\n';
  }
  std::cout << '\n' << std::endl;

  {
    std::ofstream out("harmonic-oscillator-wavefunctions.txt");
    out.precision(std::numeric_limits<data_t>::max_digits10);
    for (data_t x = static_cast<data_t>(-101) / 10;
         x <= static_cast<data_t>(101) / 10;
         x += static_cast<data_t>(1) / 100) {
      out << x << "\t" << 0.5 * x * x;
      for (const auto &r : harmonicOscillator) {
        out << "\t" << r.wavefunction(x) + r.energy;
      }
      out << "\n";
    }
  }
}

void radialHydrogen() {
  std::cout.precision(std::numeric_limits<data_t>::max_digits10);
  std::vector<Eigenspace> hydrogen = solveRadialHydrogen();

  std::cout << "Hydrogen (L=" << L << ") eigenvalues:\n\n";
  std::cout << std::setw(2) << "n"
            << "\t" << std::setw(25) << "energy"
            << "\t" << std::setw(25) << "relative deviation" << '\n';
  for (size_t i = 0; i < hydrogen.size(); i++) {
    const data_t expected =
        static_cast<data_t>(-1) / pow(static_cast<data_t>(L + i + 1), 2);
    const data_t relativeDev = abs((hydrogen[i].energy - expected) / expected);
    std::cout << std::setw(2) << i + L + 1 << "\t" << std::setw(25)
              << hydrogen[i].energy << "\t" << std::setw(25) << relativeDev
              << '\n';
  }
  std::cout << '\n' << std::endl;

  {
    std::ofstream out("radial-hydrogen-wavefunctions.txt");
    out.precision(std::numeric_limits<data_t>::max_digits10);
    for (data_t x = static_cast<data_t>(0); x <= static_cast<data_t>(1000);
         x = x * static_cast<data_t>(102) / 100 +
             static_cast<data_t>(1) / 1000) {
      out << x << "\t" << ((x < 0.75) ? -2.0 / 0.75 : -2.0 / x);
      for (const auto &r : hydrogen) {
        out << "\t" << r.wavefunction(x) + r.energy;
      }
      out << "\n";
    }
  }
}

int main(int argc, char **argv) {
  if (argc != 2 || (std::string(argv[1]) != "harmonic_oscillator" &&
                    std::string(argv[1]) != "harmonic_oscillator_spline" &&
                    std::string(argv[1]) != "hydrogen")) {
    std::cout << "Usage:\n\t" << argv[0]
              << " <example>\n\tWhere <example> is either "
                 "\"harmonic_oscillator\", \"harmonic_oscillator_spline\" or "
                 "\"hydrogen\"."
              << std::endl;
    return 1;
  }

  const std::string example(argv[1]);

  if (example == "harmonic_oscillator") {
    harmonicOscillator(false);
  } else if (example == "harmonic_oscillator_spline") {
    harmonicOscillator(true);
  } else if (example == "hydrogen") {
    radialHydrogen();
  }

  return 0;
}
