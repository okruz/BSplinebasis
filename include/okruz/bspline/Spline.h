#ifndef OKRUZ_BSPLINE_SPLINE_H
#define OKRUZ_BSPLINE_SPLINE_H
#include <okruz/bspline/exceptions/BSplineException.h>
#include <okruz/bspline/internal/misc.h>
#include <okruz/bspline/support/Support.h>

#include <algorithm>
#include <array>
#include <type_traits>
#include <vector>

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
using namespace okruz::bspline::exceptions;

/*!
 * Spline class representing spline of datatype T and order order.
 * The coefficients of the spline are defined with respect to the center point
 * xm of each interval.
 *
 * @tparam T Datatype of the spline.
 * @tparam order Order of the spline.
 */
template <typename T, size_t order>
class Spline {
 private:
  /*! Number of coefficients per interval. */
  static constexpr size_t ARRAY_SIZE = order + 1;
  /*! The support of this spline, represented by N+1 grid
                            points. */
  Support<T> _support;
  /*! Coefficients of the polynomials on each interval. */
  std::vector<std::array<T, ARRAY_SIZE>> _coefficients;

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
      return -1;  // x is not part of the spline's support

    const auto begin = _support.begin();
    const auto it = std::lower_bound(begin, _support.end(), x);

    int retval = std::distance(begin, it);
    return std::max(0, retval - 1);
  };

  /**
   * Performs consistency checks on the internal data.
   */
  void checkValidity() const {
    bool isValid =
        (!_support.containsIntervals() && _coefficients.size() == 0) ||
        (_support.size() >= 2 &&
         _coefficients.size() == _support.numberOfIntervals());
    if (!isValid) throw BSplineException(ErrorCode::INCONSISTENT_DATA);
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
   * Provides acces to the data type T of the spline.
   */
  using data_type = T;

  /*!
   * Provides access to the order of the spline.
   */
  static constexpr size_t spline_order = order;

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
  const std::vector<std::array<T, ARRAY_SIZE>> &getCoefficients()
      const noexcept {
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
    if (index < 0) return static_cast<T>(0);
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
   * empty, an exception is thrown.
   *
   * @throws BSplineException If the spline's support is empty.
   */
  const T &front() const { return _support.front(); };

  /*!
   * Returns the end of the support of this spline. If the spline is empty, an
   * exception is thrown.
   *
   * @throws BSplineException If the spline's support is empty.
   */
  const T &back() const { return _support.back(); };

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
  }

  /*!
   * Checks whether this spline returns zero for all x. Can be the case, either
   * if the support contains no intervals (i.e. the vector intervals is empty)
   * or if all coefficients are zero.
   */
  bool isZero() const {
    static const T ZERO = static_cast<T>(0);
    if (!_support.containsIntervals()) return true;
    for (const auto &cs : _coefficients) {
      for (const auto &c : cs) {
        if (c != ZERO) return false;
      }
    }
    return true;
  };

  //######################## Operator definitions ########################
  //######################################################################

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
      for (size_t j = 0; j < coeffsi.size(); j++) ncoeffsi[j] = coeffsi[j];
    }
    setData(a.getSupport(), std::move(ncoefficients));
    return *this;
  }

  /*!
   * Spline-spline multiplication operator. Returns a spline of order order +
   * ordera.
   *
   * @param a Spline to be multiplied with this spline.
   * @tparam ordera Order of spline a.
   */
  template <size_t ordera>
  Spline<T, order + ordera> operator*(const Spline<T, ordera> &a) const {
    if (!_support.hasSameGrid(a.getSupport())) {
      throw BSplineException(ErrorCode::DIFFERING_GRIDS);
    }

    static constexpr size_t NEW_ORDER = order + ordera;
    static constexpr size_t NEW_ARRAY_SIZE = NEW_ORDER + 1;

    Support newSupport = _support.calcIntersection(a.getSupport());
    const size_t nintervals = newSupport.numberOfIntervals();

    if (nintervals == 0)
      return Spline<T, NEW_ORDER>(std::move(newSupport), {});  // No overlap

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
  }

  /*!
   * Addition operator. Adds spline a to this spline.
   *
   * @param a Spline to be added.
   * @tparam ordera Order of spline a.
   */
  template <size_t ordera>
  Spline<T, std::max(order, ordera)> operator+(
      const Spline<T, ordera> &a) const {
    if (!_support.hasSameGrid(a.getSupport())) {
      throw BSplineException(ErrorCode::DIFFERING_GRIDS);
    }

    static constexpr size_t NEW_ORDER = std::max(order, ordera);
    static constexpr size_t NEW_ARRAY_SIZE = NEW_ORDER + 1;

    Support newSupport = _support.calcUnion(a.getSupport());

    const size_t nintervals = newSupport.numberOfIntervals();

    std::vector<std::array<T, NEW_ARRAY_SIZE>> ncoefficients(nintervals);

    for (size_t i = 0; i < nintervals; i++) {
      const auto absIndex = newSupport.absoluteFromRelative(i);
      const auto thisRelIndex = _support.intervalIndexFromAbsolute(absIndex);
      const auto aRelIndex = a.getSupport().intervalIndexFromAbsolute(absIndex);

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
  }

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
    (*this) = (*this) + a;
    return *this;
  }

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
  Spline<T, std::max(order, ordera)> operator-(
      const Spline<T, ordera> &a) const {
    return (*this) + (static_cast<T>(-1) * a);
  }

  //####################### Spline transformations #######################
  //######################################################################

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
        if (j > 0) coeffsi[j] += coeffs_old[j - 1];
        if (j < coeffs_old.size()) coeffsi[j] += xm * coeffs_old[j];
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
  template <size_t n = 1>
  Spline<T, orderdx(order, n)> dx() const {
    if constexpr (n > order)
      return Spline<T, 0>(Support(_support.getGrid()), {});
    else if constexpr (n == 0)
      return *this;
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
          for (size_t j = 0; j < n; j++) faculty *= i - j;
          nc[i - n] = faculty * c[i];
        }
      }
      return Spline<T, NEW_ORDER>(_support, std::move(ncoeffs));
    }
  }

  /*!
   * Calculates the second derivative.
   */
  Spline<T, orderdx(order, 2)> dx2() const { return this->template dx<2>(); };

  /*!
   * Calculates the third derivative.
   */
  Spline<T, orderdx(order, 3)> dx3() const { return this->template dx<3>(); };

};  // class Spline

//################### End of defintion of Spline class ###################
//########################################################################

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
}

/*!
 * Calculates the linear combination of splines. Is more efficient than
 * successive scalar multiplications and spline additions.
 *
 * @param coeffsBegin The iterator referencing the first element of the
 * coefficient collection.
 * @param coeffsEnd The iterator referencing the end of the coefficient
 * collection.
 * @param splinesBegin The iterator referencing the first element of the spline
 * collection.
 * @param splinesEnd The iterator referencing the end of the spline collection.
 * @tparam CeffIter An iterator referencing a coefficient of type T.
 * @tparam SplineIter An iterator referenchig a spline of type Spline<T, order>.
 * @returns The linear combination as a spline of type Spline<T, order>.
 * @throws BSplineException If the number of coefficients differs from the
 * number of splines, if the number of coefficients and splines are zero or the
 * grids of all splines are not logically equivalent.
 */
template <typename CoeffIter, typename SplineIter>
decltype(auto) linearCombination(CoeffIter coeffsBegin, CoeffIter coeffsEnd,
                                 SplineIter splinesBegin,
                                 SplineIter splinesEnd) {
  // The spline data type.
  using Spline = typename std::remove_cv_t<
      typename std::iterator_traits<SplineIter>::value_type>;

  // The coefficient data type.
  using T = typename std::remove_cv_t<
      typename std::iterator_traits<CoeffIter>::value_type>;

  // The order of the spline.
  constexpr size_t order = Spline::spline_order;

  // Check, the data type of the spline and the coefficients are consistent.
  static_assert(
      std::is_same<typename Spline::data_type, T>::value,
      "Coefficients must be of the same type as the data type of the spline.");

  {
    // The number of coefficients and splines.
    const int coeffsSize = std::distance(coeffsBegin, coeffsEnd);
    const int splinesSize = std::distance(splinesBegin, splinesEnd);

    // Check, the data is consistent.
    if (coeffsSize != splinesSize) {
      throw BSplineException(
          ErrorCode::INCONSISTENT_DATA,
          "The number of coefficients and splines must coincide.");
    }

    // std::distance() may return negative values.
    if (coeffsSize <= 0) {
      throw BSplineException(
          ErrorCode::MISSING_DATA,
          "The number of coefficients and splines may not be zero.");
    }
  }

  const auto &support0 = splinesBegin->getSupport();
  size_t startIndex = support0.getStartIndex();
  size_t endIndex = support0.getEndIndex();

  // Determine the union of all the splines' supports.
  for (auto it = splinesBegin + 1; it < splinesEnd; it++) {
    if (!it->getSupport().hasSameGrid(support0)) {
      throw BSplineException(ErrorCode::DIFFERING_GRIDS);
    }

    const size_t si = it->getSupport().getStartIndex();
    const size_t ei = it->getSupport().getEndIndex();

    if (si < startIndex) {
      startIndex = si;
    }

    if (ei > endIndex) {
      endIndex = ei;
    }
  }

  // Set up support and coefficients vector for the returned spline.
  Support newSupport(support0.getGrid(), startIndex, endIndex);
  std::vector<std::array<T, order + 1>> newCoefficients(
      newSupport.numberOfIntervals(),
      internal::make_array<T, order + 1>(static_cast<T>(0)));

  auto coeffIt = coeffsBegin;
  auto splineIt = splinesBegin;
  while (splineIt < splinesEnd) {
    // Get spline and coefficient.
    const auto &spline = *splineIt;
    const T coeff = *coeffIt;

    for (size_t j = 0; j < spline.getSupport().numberOfIntervals(); j++) {
      // Index of the interval relative to the global grid.
      const size_t absoluteIndex = spline.getSupport().absoluteFromRelative(j);

      // Index of the interval relative to the newSupport.
      const size_t newSupportIndex =
          newSupport.intervalIndexFromAbsolute(absoluteIndex).value();

      const std::array<T, order + 1> &splineCoeffs =
          spline.getCoefficients().at(j);
      std::array<T, order + 1> &newCoeffs = newCoefficients.at(newSupportIndex);

      for (size_t k = 0; k < order + 1; k++) {
        newCoeffs[k] += coeff * splineCoeffs[k];
      }
    }

    splineIt++;
    coeffIt++;
  }

  return Spline(std::move(newSupport), std::move(newCoefficients));
}

}  // namespace okruz::bspline
#endif  // OKRUZ_BSPLINE_SPLINE_H
