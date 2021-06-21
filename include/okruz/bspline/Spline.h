#ifndef OKRUZ_BSPLINE_SPLINE_H
#define OKRUZ_BSPLINE_SPLINE_H
#include <algorithm>
#include <array>
#include <assert.h>
#include <vector>

#include <okruz/bspline/internal/misc.h>
#include <okruz/bspline/support/Support.h>

/*
 * template<typename T, size_t order>
 * class Spline
 *
 * Represents a spline of dataytype T and order order. The datatype has to
 * fulfill the following requirements:
 *   - comparisons <, <=, >, >=, == and != have to be implemented.
 *   - arithmetic operators + - * /  += -= *= /= have to be implemented
 *   - A pathway must exist, such that integer values of type T can be
 * constructed via static_cast<T>(int) (e. g. via a constructor taking an int).
 * BSplines can be generated via the method generateBspline(...).
 *
 * All methods accessing two splines assume that these splines are defined on
 * the same grid (i.e. that both splines have the same interval boundaries
 * within the intersection of their respective supports). This may also cause
 * problems when splines are constructed by adding up multiple splines. To be
 * safe, make sure that the supports of two splines being added overlap at least
 * in one grid point.
 *
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

/*!
 * Main namespace for this library.
 */
namespace okruz::bspline {

using namespace support;

//####################################################### Beginning of defintion
// of Spline class
//############################################################################################
/*!
 * Spline class representing spline of datatype T and order order.
 * The coefficients of the spline are defined with respect to the center point
 * xm of each interval.
 *
 * @tparam T Datatype of the spline.
 * @tparam order Order of the spline.
 */
template <typename T, size_t order> class Spline {
private:
  static constexpr size_t ARRAY_SIZE =
      order + 1;       /*! Number of coefficients per interval. */
  Support<T> _support; /*! The support of this spline, represented by N+1 grid
                          points. */
  std::vector<std::array<T, ARRAY_SIZE>>
      _coefficients; /*! Coefficients of the polynomials on each interval. */

  /*!
   * Finds the interval in which x lies by binary search. Used during the
   * evaluation of the spline.
   *
   * @param x Point, whose interval will be searched.
   * @return The index corresponding to the beginning of the interval which
   * contains x or -1 if x is not part of the spline's support.
   */
  int findInterval(const T &x) const {
    if (_support.size() < 2 || x > _support.back() || x < _support.front())
      return -1; // x is not part of the spline's support
    int starti = 0, endi = int(_support.size()) - 1;
    while (endi - starti > 1) {
      int middlei = (endi + starti) / 2;
      if (x > _support[middlei])
        starti = middlei;
      else
        endi = middlei;
    }
    return starti;
  };

  /**
   * Performs consistency checks on the internal data.
   */
  void checkValidity() const {
    assert((!_support.containsIntervals() && _coefficients.size() == 0) ||
           (_support.size() >= 2 &&
            _coefficients.size() == _support.numberOfIntervals()));
  };

  /*!
   * Resets the data of the spline and performs sanity checks.
   *
   * @param intervals The grid points representing the intervals on which the
   * spline is defined.
   * @param coefficients Polynomial coefficients on each of the intervals.
   */
  void setData(Support<T> support,
               std::vector<std::array<T, ARRAY_SIZE>> coefficients) {
    _support = std::move(support);
    _coefficients = std::move(coefficients);
    checkValidity();
  };

public:
  /*!
   * Constructor setting the data. Performs sanity checks.
   *
   * @param support The spline's support.
   * @param coefficients Polynomial coefficients of the spline on each interval.
   */
  Spline(Support<T> support,
         std::vector<std::array<T, ARRAY_SIZE>> coefficients)
      : _support(std::move(support)), _coefficients(std::move(coefficients)) {
    checkValidity();
  };

  /**
   * Constructs an empty spline on the global grid.
   *
   * @param grid The global grid.
   */
  Spline(Grid<T> grid) : Spline(Support(std::move(grid)), {}){};

  // Default constructors and operators generated by the compiler.
  Spline() = delete;
  Spline(const Spline &m) = default;
  Spline(Spline &&m) = default;
  virtual ~Spline() = default;
  Spline &operator=(const Spline &m) = default;
  Spline &operator=(Spline &&m) = default;

  /*!
   * Returns the spline's support.
   */
  const Support<T> &getSupport() const noexcept { return _support; };

  /*!
   * Returns the polynomial coefficients of the spline for each interval.
   */
  const std::vector<std::array<T, ARRAY_SIZE>> &
  getCoefficients() const noexcept {
    return _coefficients;
  };

  /*!
   * Evaluates the spline at point x.
   *
   * @param x Point at which to evaluate the spline. If x is outside of the
   * support of the spline, zero is returned.
   */
  T operator()(const T &x) const {
    int index = findInterval(x);
    if (index < 0)
      return static_cast<T>(0);
    const auto &coeffs = _coefficients[index];

    // distance between x and the middlepoint of the interval
    const T dx =
        x - (_support[index + 1] + _support[index]) / static_cast<T>(2);
    T xpot = static_cast<T>(1), result = static_cast<T>(0);
    for (const T &c : coeffs) {
      result += xpot * c;
      xpot *= dx;
    }
    return result;
  };

  /*!
   * Returns the beginning of the support of this spline. If the spline is
   * empty, zero is returned.
   */
  T start() const {
    if (_support.size() == 0) {
      return static_cast<T>(0);
    }
    return _support.front();
  };

  /*!
   * Returns the end of the support of this spline. If the spline is empty, zero
   * is returned.
   */
  T end() const {
    if (_support.size() == 0) {
      return static_cast<T>(0);
    }
    return _support.back();
  };

  /*!
   * Checks whether the supports of the two splines overlap.
   *
   * @param m2 Other spline against which to check.
   * @tparam order2 Order of spline m2.
   */
  template <size_t order2>
  bool checkOverlap(const Spline<T, order2> &m2) const {
    if (!_support.containsIntervals() || !m2.getSupport().containsIntervals())
      return false;
    const bool isNotOverlapping = m2.getSupport().back() <= _support.front() ||
                                  m2.getSupport().front() >= _support.back();
    return !isNotOverlapping;
  };

  /*!
   * Checks whether this spline returns zero for all x. Can be the case, either
   * if the support contains no intervals (i.e. the vector intervals is empty)
   * or if all coefficients are zero.
   */
  bool isZero() const {
    if (!_support.containsIntervals())
      return true;
    const T zero = static_cast<T>(0);
    for (const auto &cs : _coefficients) {
      for (const auto &c : cs) {
        if (c != zero)
          return false;
      }
    }
    return true;
  };

  // ################################ Operator definitions
  // ###############################################

  /*!
   * Scalar-division operator. Divides this spline by the scalar d.
   *
   * @param d Scalar by which to divide this spline.
   */
  Spline<T, order> operator/(const T &d) const {
    return (*this) * (static_cast<T>(1) / d);
  };

  /*!
   * Scalar-multiplication operator. Multiplies this spline with the scalar d.
   *
   * @param d Scalar by which to multiply this spline.
   */
  Spline<T, order> operator*(const T &d) const {
    Spline<T, order> ret(*this);
    for (auto &cs : ret._coefficients) {
      for (auto &c : cs) {
        c *= d;
      }
    }
    return ret;
  };

  /*!
   * In-place scalar-multiplication operator. Multiplies this spline with the
   * scalar d in-place.
   *
   * @param d Scalar by which to multiply this spline.
   */
  Spline<T, order> &operator*=(const T &d) {
    for (auto &cs : _coefficients) {
      for (auto &c : cs) {
        c *= d;
      }
    }
    return *this;
  };

  /*!
   * In-place scalar-division operator. Divides this spline by the scalar d
   * in-place.
   *
   * @param d Scalar by which to divide this spline.
   */
  Spline<T, order> &operator/=(const T &d) {
    (*this) *= (static_cast<T>(1) / d);
    return *this;
  };

  /*!
   * Unary minus operator.
   */
  Spline<T, order> operator-() const { return (*this) * static_cast<T>(-1); };

  /*!
   * Copy assign of spline to this spline object. The operation is only well
   * defined if the order of the spline to be assigned is lower than or equal to
   * the order of this spline object.
   *
   * @param a Spline to be assigned.
   * @tparam ordera Order of spline a.
   */
  template <size_t ordera>
  Spline<T, order> &operator=(const Spline<T, ordera> &a) {
    // The case ordera == order should be handled by the default assignment
    // operator which is automatically generated.
    static_assert(
        ordera < order,
        "The assignment operator is only defined if the order of the rhs "
        "spline is lower than or equal to that of the lhs spline.");

    std::vector<std::array<T, ARRAY_SIZE>> ncoefficients(
        a.getCoefficients().size(),
        internal::make_array<T, ARRAY_SIZE>(static_cast<T>(0)));
    for (size_t i = 0; i < a.getCoefficients().size(); i++) {
      const auto &coeffsi = a.getCoefficients()[i];
      auto &ncoeffsi = ncoefficients[i];
      for (size_t j = 0; j < coeffsi.size(); j++)
        ncoeffsi[j] = coeffsi[j];
    }
    setData(a.getSupport(), std::move(ncoefficients));
    return *this;
  };

  /*!
   * Spline-spline multiplication operator. Returns a spline of order order +
   * ordera.
   *
   * @param a Spline to be multiplied with this spline.
   * @tparam ordera Order of spline a.
   */
  template <size_t ordera>
  Spline<T, order + ordera> operator*(const Spline<T, ordera> &a) const {
    assert(_support.hasSameGrid(a.getSupport()));

    static constexpr size_t NEW_ORDER = order + ordera;
    static constexpr size_t NEW_ARRAY_SIZE = NEW_ORDER + 1;

    Support newSupport = _support.calcIntersection(a.getSupport());
    const size_t nintervals = newSupport.numberOfIntervals();

    if (nintervals == 0)
      return Spline<T, NEW_ORDER>(std::move(newSupport), {}); // No overlap

    std::vector<std::array<T, NEW_ARRAY_SIZE>> newCoefficients(
        nintervals, internal::make_array<T, NEW_ARRAY_SIZE>(static_cast<T>(0)));

    for (size_t i = 0; i < nintervals; i++) {
      const auto ai = newSupport.absoluteFromRelative(i);

      auto thisRelIndex = _support.relativeFromAbsolute(ai).value();
      auto aRelIndex = a.getSupport().relativeFromAbsolute(ai).value();

      const auto &thiscoeffs = _coefficients[thisRelIndex];
      const auto &acoeffs = a.getCoefficients()[aRelIndex];
      auto &coeffsi = newCoefficients[i];

      for (size_t j = 0; j < order + 1; j++) {
        for (size_t k = 0; k < ordera + 1; k++) {
          coeffsi[j + k] += thiscoeffs[j] * acoeffs[k];
        }
      }
    }
    return Spline<T, NEW_ORDER>(std::move(newSupport),
                                std::move(newCoefficients));
  };

  /*!
   * Addition operator. Adds spline a to this spline.
   *
   * @param a Spline to be added.
   * @tparam ordera Order of spline a.
   */
  template <size_t ordera>
  Spline<T, std::max(order, ordera)>
  operator+(const Spline<T, ordera> &a) const {
    assert(_support.hasSameGrid(a.getSupport()));

    static constexpr size_t NEW_ORDER = std::max(order, ordera);
    static constexpr size_t NEW_ARRAY_SIZE = NEW_ORDER + 1;

    Support newSupport = _support.calcUnion(a.getSupport());

    const size_t nintervals = newSupport.numberOfIntervals();

    std::vector<std::array<T, NEW_ARRAY_SIZE>> ncoefficients(nintervals);

    for (size_t i = 0; i < nintervals; i++) {
      const auto absIndex = newSupport.absoluteFromRelative(i);
      const auto thisRelIndex = _support.relativeFromAbsolute(absIndex);
      const auto aRelIndex = a.getSupport().relativeFromAbsolute(absIndex);

      if (thisRelIndex && !aRelIndex) {
        ncoefficients[i] =
            internal::changearraysize<T, order + 1, NEW_ARRAY_SIZE>(
                _coefficients[thisRelIndex.value()]);
      } else if (aRelIndex && !thisRelIndex) {
        ncoefficients[i] =
            internal::changearraysize<T, ordera + 1, NEW_ARRAY_SIZE>(
                a.getCoefficients()[aRelIndex.value()]);
      } else if (thisRelIndex && aRelIndex) {
        ncoefficients[i] = internal::add<T, ordera + 1, order + 1>(
            a.getCoefficients()[aRelIndex.value()],
            _coefficients[thisRelIndex.value()]);
      } else {
        ncoefficients[i] =
            internal::make_array<T, NEW_ARRAY_SIZE>(static_cast<T>(0));
      }
    }
    return Spline<T, NEW_ORDER>(std::move(newSupport),
                                std::move(ncoefficients));
  };

  /*!
   * In-place addition operator. Adds spline a to this spline. The operation is
   * only well defined if the order of spline a is lower than or equal to the
   * order of this spline object.
   *
   * @param a Spline to be added.
   * @tparam ordera Order of spline a.
   */
  template <size_t ordera>
  Spline<T, order> &operator+=(const Spline<T, ordera> &a) {
    static_assert(
        ordera <= order,
        "The operators += and -= are only defined if the order of the rhs "
        "spline is lower than or equal to that of the lhs spline.");
    assert(_support.hasSameGrid(a.getSupport()));

    static constexpr size_t NEW_ORDER = std::max(order, ordera);
    static constexpr size_t NEW_ARRAY_SIZE = NEW_ORDER + 1;

    Support newSupport = _support.calcUnion(a.getSupport());

    const size_t nintervals = newSupport.numberOfIntervals();

    std::vector<std::array<T, NEW_ARRAY_SIZE>> ncoefficients(nintervals);

    for (size_t i = 0; i < nintervals; i++) {
      const auto absIndex = newSupport.absoluteFromRelative(i);
      const auto thisRelIndex = _support.relativeFromAbsolute(absIndex);
      const auto aRelIndex = a.getSupport().relativeFromAbsolute(absIndex);

      if (thisRelIndex && !aRelIndex) {
        ncoefficients[i] =
            internal::changearraysize<T, order + 1, NEW_ARRAY_SIZE>(
                _coefficients[thisRelIndex.value()]);
      } else if (aRelIndex && !thisRelIndex) {
        ncoefficients[i] =
            internal::changearraysize<T, ordera + 1, NEW_ARRAY_SIZE>(
                a.getCoefficients()[aRelIndex.value()]);
      } else if (thisRelIndex && aRelIndex) {
        ncoefficients[i] = internal::add<T, ordera + 1, order + 1>(
            a.getCoefficients()[aRelIndex.value()],
            _coefficients[thisRelIndex.value()]);
      } else {
        ncoefficients[i] =
            internal::make_array<T, NEW_ARRAY_SIZE>(static_cast<T>(0));
      }
    }
    setData(std::move(newSupport), std::move(ncoefficients));
    return *this;
  };

  /*!
   * Binary in-place subtraction operator. Subtracts spline a from this spline.
   * The operation is only well defined if the order of the spline to be
   * subtracted is lower than or equal to the order of this spline object.
   *
   * @param a Spline to be subtracted.
   * @tparam ordera Order of spline a.
   */
  template <size_t ordera>
  Spline<T, order> &operator-=(const Spline<T, ordera> &a) {
    (*this) += (static_cast<T>(-1) * a);
    return *this;
  }

  /*!
   * Binary subtraction operator. Subtracts spline a from this spline.
   *
   * @param a Spline to be subtracted.
   * @tparam ordera Order of spline a.
   */
  template <size_t ordera>
  Spline<T, std::max(order, ordera)>
  operator-(const Spline<T, ordera> &a) const {
    return (*this) + (static_cast<T>(-1) * a);
  }

  // ################################### Spline transformations
  // ###########################################################

  /*!
   * Returns a spline g(x) = x f(x), where f(x) is this spline.
   */
  Spline<T, order + 1> timesx() const {
    std::vector<std::array<T, ARRAY_SIZE + 1>> newcoeffs(
        _coefficients.size(),
        internal::make_array<T, ARRAY_SIZE + 1>(static_cast<T>(0)));
    for (size_t i = 0; i + 1 < _support.size(); i++) {
      const T xm = (_support[i + 1] + _support[i]) / static_cast<T>(2);
      const std::array<T, ARRAY_SIZE> &coeffs_old = _coefficients[i];
      auto &coeffsi = newcoeffs[i];
      for (size_t j = 0; j <= coeffs_old.size(); j++) {
        if (j > 0)
          coeffsi[j] += coeffs_old[j - 1];
        if (j < coeffs_old.size())
          coeffsi[j] += xm * coeffs_old[j];
        ;
      }
    }
    return Spline<T, order + 1>(_support, std::move(newcoeffs));
  };

  /*!
   * Calculates the order of a spline representing a derivative of another
   * spline. Defined for convenience.
   *
   * @param spline_order Order of the spline before applying the derivative.
   * @param derivative_order Order of the derivative to be applied.
   */
  static constexpr size_t orderdx(size_t spline_order,
                                  size_t derivative_order) {
    if (derivative_order > spline_order)
      return 0;
    else
      return spline_order - derivative_order;
  }

  /*!
   * Returns a spline g(x) = \\frac{\\partial^n}{\\partial x^n} f(x), where f(x)
   * is this spline. Assumes the spline is n-1 times continously differentiable.
   *
   * @tparam n Order of the derivative.
   */
  template <size_t n = 1> Spline<T, orderdx(order, n)> dx() const {
    if constexpr (n > order)
      return Spline<T, 0>(Support(_support.getGrid()), {});
    else {
      static constexpr size_t NEW_ORDER = orderdx(order, n);
      static constexpr size_t NEW_ARRAY_SIZE = NEW_ORDER + 1;
      std::vector<std::array<T, NEW_ARRAY_SIZE>> ncoeffs(
          _coefficients.size(),
          internal::make_array<T, NEW_ARRAY_SIZE>(static_cast<T>(0)));
      for (size_t ii = 0; ii < _coefficients.size(); ii++) {
        auto &nc = ncoeffs[ii];
        const auto &c = _coefficients[ii];
        for (size_t i = n; i < c.size(); i++) {
          size_t faculty = 1;
          for (size_t j = 0; j < n; j++)
            faculty *= i - j;
          nc[i - n] = faculty * c[i];
        }
      }
      return Spline<T, NEW_ORDER>(_support, std::move(ncoeffs));
    }
  };

  /*!
   * Calculates the second derivative.
   */
  Spline<T, orderdx(order, 2)> dx2() const { return this->template dx<2>(); };

  /*!
   * Calculates the third derivative.
   */
  Spline<T, orderdx(order, 3)> dx3() const { return this->template dx<3>(); };

}; // class Spline
//####################################################### End of defintion of
// Spline class
//############################################################################################

/*!
 * Commutation of spline scalar multiplication operator.
 *
 * @param d Scalar to be multiplied.
 * @param b Spline to be multiplied.
 * @tparam T Datatype of spline and scalar.
 * @tparam order Order of the spline.
 */
template <typename T, size_t order>
inline Spline<T, order> operator*(const T &d, const Spline<T, order> &b) {
  return b * d;
};

};     // namespace okruz::bspline
#endif // OKRUZ_BSPLINE_SPLINE_H
