#ifndef MYSPLINE_H
#define MYSPLINE_H
#include <vector>
#include <assert.h>
#include <algorithm>
#include <array>

/*
 * template<typename T, size_t order>
 * class myspline::myspline
 * 
 * Represents a spline of dataytype T and order order. The datatype has to fulfill the following requirements:
 *   - comparisons <, <=, >, >=, == and != have to be implemented.
 *   - arithmetic operators + - * /  += -= *= /= have to be implemented
 *   - A pathway must exist, such that integer values of type T can be constructed via static_cast<T>(int) (e. g. via a constructor taking an int).
 * BSplines can be generated via the method generateBspline(...).
 * 
 * All methods accessing two splines assume that these splines are defined on the same grid (i.e. that
 * both splines have the same interval boundaries within the intersection of their respective supports).
 * This may also cause problems when splines are constructed by adding up
 * multiple splines. To be safe, make sure that the supports of two splines being added overlap at least in one grid point.
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
namespace myspline {

template<typename T, size_t order>
class myspline ;

namespace internal {
template<typename T, size_t order1, size_t order2>
void findOverlappingIntervals(const myspline<T, order1> &m1, const myspline<T, order2> &m2, size_t &startindex1, size_t &startindex2, size_t &nintervals);

/*!
 * Creates an std::array<T, size> with all values set to val.
 *
 * @param val Value to initialise the members of the array with.
 * @tparam T Datatype of the array.
 * @tparam size Size of the array.
 */
template<typename T, size_t size>
std::array<T, size> make_array(T val) {
    std::array<T, size> ret;
    ret.fill(val);
    return ret;
}

/*!
 * Adds two vectors, taking into account that one vector might be longer than the other one.
 *
 * @param a First array to be added.
 * @param b Second array to be added.
 * @tparam T Datatype of both arrays.
 * @tparam sizea Size of the first array.
 * @tparam sizeb Size of the second array.
 */
template<typename T, size_t sizea, size_t sizeb>
std::array<T, std::max(sizea, sizeb)> add(const std::array<T, sizea> &a, const std::array<T, sizeb> &b) {
    static constexpr size_t NEW_ARRAY_SIZE = std::max(sizea, sizeb);
    if constexpr (sizeb > sizea){
        return add(b, a);
    } else {
        std::array<T, NEW_ARRAY_SIZE> ret = a;
        for(size_t i = 0; i < sizeb; i++) {
            ret[i] += b[i];
        }
        return ret;
    }
};

/*!
 * Copies one array into a larger one, filling up the additional members with zeros.
 * 
 * @param in Array to be copied.
 * @tparam T Datatype of the arrays.
 * @tparam sizein Size of the input array.
 * @tparam sizeout Size of the output array. Must fulfil sizeout >= sizein.
 */
template<typename T, size_t sizein, size_t sizeout>
std::array<T, sizeout> changearraysize(const std::array<T, sizein> &in) {
    static_assert(sizeout >= sizein, "sizeout must be bigger or equal to sizein.");
    if constexpr(sizeout == sizein) return in;
    else {
        std::array<T, sizeout> ret = make_array<T, sizeout>(static_cast<T>(0));
        for (size_t i = 0; i < sizein; i++) ret[i] = in[i];
        return ret;
    }
}

/*!
 * Removes double entries from the vector list and sorts it.
 *
 * @param list The list to be manipulated.
 * @tparam T Datatype of the list.
 */
template<typename T>
void makeuniquesorted(std::vector<T> &list) {
    std::sort(list.begin(), list.end());
    auto ip = std::unique(list.begin(), list.end());
    list = std::vector<T>(list.begin(), ip);
};


/*!
 * Checks whether the elements of the vector are are steadily increasing.
 * 
 * @param vec Vector to be checked.
 * @tparam T Datatype of the elements of the vector vec.
 */
template <typename T>
bool isSteadilyIncreasing(const std::vector<T> &vec) {
    for(size_t i = 1; i < vec.size(); i++) {
        if (vec[i-1] >= vec[i]) return false;
    }
    return true;
};
}; // end namespace internal
//####################################################### Beginning of defintion of myspline class ############################################################################################
/*!
 * Spline class representing spline of datatype T and order order.
 * The coefficients of the spline are defined with respect to the center point xm of each interval.
 *
 * @tparam T Datatype of the spline.
 * @tparam order Order of the spline.
 */
template<typename T, size_t order>
class myspline {
    private:
        static constexpr size_t ARRAY_SIZE = order+1;                 /*! Number of coefficients per interval. */
        std::vector<T> _intervals;                                    /*! The N intervals of the support of this spline, represented by N+1 grid points. */
        std::vector<std::array<T, ARRAY_SIZE>> _coefficients;         /*! Coefficients of the polynomials on each interval. */


        /*!
         * Finds the interval in which x lies by binary search. Used during the evaluation of the spline.
         *
         * @param x Point, whose interval will be searched.
         * @return The index corresponding to the beginning of the interval which contains x or -1 if x is not part of the spline's support.
         */
        int findInterval(const T& x) const {
            if(_intervals.size() < 2 || x > _intervals.back() || x < _intervals.front()) return -1; // x is not part of the spline's support
            int starti = 0, endi = int(_intervals.size())-1;
            while(endi-starti > 1) {
                int middlei = (endi+starti)/2;
                if (x > _intervals[middlei]) starti = middlei;
                else endi = middlei;
            }
            return starti;
        };


       /*!
        * Resets the data of the spline and performs sanity checks.
        *
        * @param intervals The grid points representing the intervals on which the spline is defined.
        * @param coefficients Polynomial coefficients on each of the intervals.
        */
       void setData(std::vector<T> intervals, std::vector<std::array<T, ARRAY_SIZE>> coefficients) {
           _intervals = std::move(intervals);
           _coefficients = std::move(coefficients);
           assert((_intervals.size() == 0 && _coefficients.size() == 0) || (_intervals.size() >= 2 && _coefficients.size() + 1 == _intervals.size()));
           assert(internal::isSteadilyIncreasing(_intervals));   
       };

   public:
       /*!
        * Constructor setting the data. Performs sanity checks.
        *
        * @param intervals Grid points symbolising the intervals of the spline's support.
        * @param coefficients Polynomial coefficients of the spline on each interval.
        */
       myspline(std::vector<T> intervals, std::vector<std::array<T, ARRAY_SIZE>> coefficients):  _intervals(std::move(intervals)), _coefficients(std::move(coefficients)) {
            assert((_intervals.size() == 0 && _coefficients.size() == 0) || (_intervals.size() >= 2 && _coefficients.size() + 1 == _intervals.size()));
            assert(internal::isSteadilyIncreasing(_intervals));    
       };

       // Default constructors and operators generated by the compiler.
       myspline() = default;
       myspline(const myspline &m) = default;
       myspline(myspline &&m) = default;
       virtual ~myspline() = default;
       myspline& operator=(const myspline &m) = default;
       myspline& operator=(myspline &&m) = default;

       /*!
        * Returns the vector of the grid points representing the intervals of the spline's support.
        */
       const std::vector<T>& getIntervals() const noexcept {return _intervals;};

       /*!
        * Returns the polynomial coefficients of the spline for each interval.
        */
       const std::vector<std::array<T,ARRAY_SIZE>>& getCoefficients() const noexcept {return _coefficients;};


      /*!
       * Evaluates the spline at point x.
       *
       * @param x Point at which to evaluate the spline. If x is outside of the support of the spline, zero is returned.
       */
       T operator()(const T &x) const {
           int index = findInterval(x);
           if (index < 0) return static_cast<T>(0);
           const auto &coeffs = _coefficients[index];

           // distance between x and the middlepoint of the interval
           const T dx = x - (_intervals[index+1] + _intervals[index])/static_cast<T>(2);
           T xpot = static_cast<T>(1), result = static_cast<T>(0);
           for (const T& c: coeffs) {
               result += xpot * c;
               xpot *= dx;
           }
           return result;
       };


        /*!
         * Returns the beginning of the support of this spline. If the spline is empty, zero is returned.
         */
        T start() const {
            if (_intervals.size() == 0) {
                return static_cast<T>(0);
            }
            return _intervals.front();
        };

        /*!
         * Returns the end of the support of this spline. If the spline is empty, zero is returned.
         */
        T end() const {
            if (_intervals.size() == 0) {
                return static_cast<T>(0);
            }
            return _intervals.back();
        };


      /*!
       * Checks whether the supports of the two splines overlap.
       *
       * @param m2 Other spline against which to check.
       * @tparam order2 Order of spline m2.
       */
       template<size_t order2>
       bool checkOverlap(const myspline<T, order2> &m2) const {
            if (_intervals.empty() || m2.getIntervals().empty()) return false;
            const bool isNotOverlapping = m2.getIntervals().back() <= _intervals.front() || m2.getIntervals().front() >= _intervals.back();
            return ! isNotOverlapping;
       };

      /*!
       * Checks whether this spline returns zero for all x. Can be the case, either if the support contains no intervals (i.e. the vector intervals is empty) or if all coefficients are zero.
       */
       bool isZero() const {
           if(_intervals.size() == 0) return true;
           const T zero = static_cast<T>(0);
           for (const auto &cs: _coefficients) {
               for(const auto& c: cs) {
                   if (c != zero) return false;
               }
           }
           return true;
       };



      // ################################ Operator definitions ###############################################

      /*!
       * Scalar-division operator. Divides this spline by the scalar d.
       *
       * @param d Scalar by which to divide this spline.
       */
       myspline<T, order> operator/(const T& d) const {
           return (*this) * (static_cast<T>(1)/d);
       };

      /*!
       * Scalar-multiplication operator. Multiplies this spline with the scalar d.
       *
       * @param d Scalar by which to multiply this spline.
       */
       myspline<T, order> operator*(const T& d) const {
            myspline<T, order> ret(*this);
            for(auto &cs: ret._coefficients){
                for (auto &c: cs) {
                    c *= d;
                }
            }
            return ret;
       };

      /*!
       * In-place scalar-multiplication operator. Multiplies this spline with the scalar d in-place.
       *
       * @param d Scalar by which to multiply this spline.
       */
       myspline<T, order>& operator*=(const T& d) {
            for(auto &cs: _coefficients){
                for (auto &c: cs) {
                    c *= d;
                }
            }
            return *this;
       };

      /*!
       * In-place scalar-division operator. Divides this spline by the scalar d in-place.
       *
       * @param d Scalar by which to divide this spline.
       */
       myspline<T, order>& operator/=(const T& d) {
           (*this) *= (static_cast<T>(1)/d);
           return *this;
       };


      /*!
       * Unary minus operator.
       */
       myspline<T, order> operator-() const {
            return (*this) * static_cast<T>(-1);
       };



      /*!
       * Copy assign of spline to this spline object. The operation is only well defined if the order of the spline to be assigned is lower than or equal to the order of this spline object.
       *
       * @param a Spline to be assigned.
       * @tparam ordera Order of spline a.
       */
       template<size_t ordera>
       myspline<T, order>& operator=(const myspline<T, ordera> &a) {
           // The case ordera == order should be handled by the default assignment operator which is automatically generated.
           static_assert(ordera < order, "The assignment operator is only defined if the order of the rhs spline is lower than or equal to that of the lhs spline.");
        
           std::vector<std::array<T, ARRAY_SIZE>> ncoefficients(a.getCoefficients().size(),  internal::make_array<T,ARRAY_SIZE>(static_cast<T>(0)));
           for (size_t i = 0; i  < a.getCoefficients().size(); i++) {
               const auto &coeffsi = a.getCoefficients()[i];
               auto &ncoeffsi = ncoefficients[i];
               for(size_t j = 0; j < coeffsi.size(); j++) ncoeffsi[j] = coeffsi[j];
           }
           setData(a.getIntervals(), std::move(ncoefficients));
           return *this;
       };


      /*!
       * Spline-spline multiplication operator. Returns a spline of order order + ordera.
       *
       * @param a Spline to be multiplied with this spline.
       * @tparam ordera Order of spline a.
       */
       template<size_t ordera>
       myspline<T, order+ ordera> operator*(const myspline<T, ordera> &a) const {
           static constexpr size_t NEW_ORDER =order + ordera;
           static constexpr size_t NEW_ARRAY_SIZE = NEW_ORDER + 1;
           size_t startindex1, startindex2, nintervals;
           internal::findOverlappingIntervals(*this, a, startindex1, startindex2, nintervals);
           if (nintervals == 0) return myspline<T, NEW_ORDER>(); // No overlap

           std::vector<std::array<T,NEW_ARRAY_SIZE>> new_coefficients(nintervals, internal::make_array<T,NEW_ARRAY_SIZE>(static_cast<T>(0)));
           std::vector<T> new_intervals(nintervals+1);

           new_intervals[0]= _intervals[startindex1];

           for (size_t i = 0; i < nintervals; i++) {
               const auto &thiscoeffs = _coefficients[i + startindex1];
               const auto &acoeffs = a.getCoefficients()[i + startindex2];
               auto &coeffsi = new_coefficients[i];
               
               new_intervals[i+1]= _intervals[startindex1 + i + 1];

               for (size_t j= 0; j < order+1; j++) {
                   for (size_t k = 0; k < ordera+1; k++) {
                       coeffsi[j+k] += thiscoeffs[j] * acoeffs[k];
                   }
               }
           }
           return myspline<T, NEW_ORDER>(std::move(new_intervals), std::move(new_coefficients));
       };

      /*!
       * Addition operator. Adds spline a to this spline. 
       *
       * @param a Spline to be added.
       * @tparam ordera Order of spline a.
       */
       template<size_t ordera>
       myspline<T, std::max(order, ordera)> operator+(const myspline<T, ordera> &a) const {
           static constexpr size_t  NEW_ORDER = std::max(order, ordera);
           static constexpr size_t NEW_ARRAY_SIZE = NEW_ORDER + 1;
           std::vector<T> nintervals; nintervals.reserve(a.getIntervals().size() + _intervals.size());
           nintervals.insert(nintervals.end(), a.getIntervals().begin(), a.getIntervals().end());
           nintervals.insert(nintervals.end(), _intervals.begin(), _intervals.end());
           internal::makeuniquesorted(nintervals);
        
           std::vector<std::array<T, NEW_ARRAY_SIZE>> ncoefficients;
           if (nintervals.size() > 1) ncoefficients.reserve(nintervals.size() - 1);

           for (size_t i = 0; i + 1 < nintervals.size(); i++) {
               const size_t posthis = std::distance(_intervals.begin(), std::find(_intervals.begin(), _intervals.end(), nintervals[i]));
               const size_t posa = std::distance(a.getIntervals().begin(), std::find(a.getIntervals().begin(), a.getIntervals().end(), nintervals[i]));
            
               const bool thisexists = (posthis < _coefficients.size());
               const bool aexists = (posa < a.getCoefficients().size());

               if (thisexists && ! aexists) {
                   ncoefficients.push_back(internal::changearraysize<T ,order +1, NEW_ARRAY_SIZE>(_coefficients[posthis]));
               } else if (aexists && ! thisexists) {
                   ncoefficients.push_back(internal::changearraysize<T, ordera + 1, NEW_ARRAY_SIZE>(a.getCoefficients()[posa]));
               } else if (thisexists && aexists){
                   ncoefficients.push_back(internal::add<T, ordera + 1, order + 1>(a.getCoefficients()[posa], _coefficients[posthis]));
               } else {
                   ncoefficients.push_back(internal::make_array<T, NEW_ARRAY_SIZE>(static_cast<T>(0)));
               }
           }
           return myspline<T, NEW_ORDER> (std::move(nintervals), std::move(ncoefficients));
       };

      /*!
       * In-place addition operator. Adds spline a to this spline. The operation is only well defined if the order of spline a is lower than or equal to the order of this spline object.
       *
       * @param a Spline to be added.
       * @tparam ordera Order of spline a.
       */
       template<size_t ordera>
       myspline<T, order>& operator+=(const myspline<T, ordera> &a) {
           static_assert(ordera <= order, "The operators += and -= are only defined if the order of the rhs spline is lower than or equal to that of the lhs spline.");
           std::vector<T> nintervals; nintervals.reserve(a.getIntervals().size() + _intervals.size());
           nintervals.insert(nintervals.end(), a.getIntervals().begin(), a.getIntervals().end());
           nintervals.insert(nintervals.end(), _intervals.begin(), _intervals.end());
           internal::makeuniquesorted(nintervals);
        
           std::vector<std::array<T, ARRAY_SIZE>> ncoefficients;
           if (nintervals.size() > 1) ncoefficients.reserve(nintervals.size() - 1);

           for (size_t i = 0; i + 1 < nintervals.size(); i++) {
               const size_t posthis = std::distance(_intervals.begin(), std::find(_intervals.begin(), _intervals.end(), nintervals[i]));
               const size_t posa = std::distance(a.getIntervals().begin(), std::find(a.getIntervals().begin(), a.getIntervals().end(), nintervals[i]));
            
               const bool thisexists = (posthis < _coefficients.size());
               const bool aexists = (posa < a.getCoefficients().size());

               if (thisexists && ! aexists) {
                   ncoefficients.push_back(_coefficients[posthis]);
               } else if (aexists && ! thisexists) {
                   ncoefficients.push_back(internal::changearraysize<T, ordera + 1, ARRAY_SIZE>(a.getCoefficients()[posa]));
               } else if (thisexists && aexists){
                   ncoefficients.push_back(internal::add<T, ordera + 1, ARRAY_SIZE>(a.getCoefficients()[posa], _coefficients[posthis]));
               } else {
                   ncoefficients.push_back(internal::make_array<T, ARRAY_SIZE>(static_cast<T>(0)));
               }
           }
           
           setData(std::move(nintervals), std::move(ncoefficients));
           return *this;
       };

      /*!
       * Binary in-place subtraction operator. Subtracts spline a from this spline. The operation is only well defined if the order of the spline to be subtracted is lower than or equal to the order of this spline object.
       *
       * @param a Spline to be subtracted.
       * @tparam ordera Order of spline a.
       */
       template<size_t ordera>
       myspline<T, order>& operator-=(const myspline<T, ordera> &a) {
           (*this) += (static_cast<T>(-1) * a);
           return *this;
       }


      /*!
       * Binary subtraction operator. Subtracts spline a from this spline. 
       *
       * @param a Spline to be subtracted.
       * @tparam ordera Order of spline a.
       */
       template<size_t ordera>
       myspline<T, std::max(order, ordera)> operator-(const myspline<T, ordera> &a) const {
           return (*this) + (static_cast<T>(-1) * a);
       }

      // ################################### Spline transformations ###########################################################


      /*!
       * Returns a spline g(x) = x f(x), where f(x) is this spline.
       */
       myspline<T, order + 1> timesx() const {
           std::vector<std::array<T, ARRAY_SIZE + 1>> newcoeffs(_coefficients.size(), internal::make_array<T, ARRAY_SIZE + 1>(static_cast<T>(0)));
           for (size_t i = 0; i+1 < _intervals.size(); i++) {
               const T xm = (_intervals[i+1] + _intervals[i])/static_cast<T>(2);
               const std::array<T, ARRAY_SIZE> &coeffs_old = _coefficients[i];
               auto &coeffsi = newcoeffs[i];
               for(size_t j = 0; j <= coeffs_old.size(); j++) {
                   if (j > 0) coeffsi[j] += coeffs_old[j-1];
                   if (j < coeffs_old.size()) coeffsi[j] += xm * coeffs_old[j];;
               }
           }
           return myspline<T, order + 1>(_intervals, std::move(newcoeffs));
       };


      /*!
       * Calculates the order of a spline representing a derivative of another spline. Defined for convenience.
       *
       * @param spline_order Order of the spline before applying the derivative.
       * @param derivative_order Order of the derivative to be applied.
       */
      static constexpr size_t orderdx(size_t spline_order, size_t derivative_order){
          if (derivative_order > spline_order) return 0;
          else return spline_order -derivative_order;
      }

      /*!
       * Returns a spline g(x) = \\frac{\\partial^n}{\\partial x^n} f(x), where f(x) is this spline.
       * Assumes the spline is n-1 times continously differentiable.
       *
       * @tparam n Order of the derivative.
       */
       template<size_t n=1>
       myspline<T, orderdx(order, n)> dx() const {
           if constexpr(n > order) return myspline<T, 0>();
           else {
               static constexpr size_t  NEW_ORDER = orderdx(order, n);
               static constexpr size_t NEW_ARRAY_SIZE = NEW_ORDER + 1;
               std::vector<std::array<T, NEW_ARRAY_SIZE>> ncoeffs(_coefficients.size(), internal::make_array<T, NEW_ARRAY_SIZE>(static_cast<T>(0)));
               for (size_t ii = 0; ii < _coefficients.size(); ii++) {
                   auto &nc = ncoeffs[ii];
                   const auto &c = _coefficients[ii];
                   for (size_t i = n; i < c.size(); i++) {
                       size_t faculty = 1;
                       for(size_t j = 0; j < n; j++) faculty *= i-j;
                       nc[i-n] = faculty * c[i];
                   }
               }
               return myspline<T, NEW_ORDER>(_intervals, std::move(ncoeffs));
           }
       };

       /*!
        * Calculates the second derivative.
        */
       myspline<T, orderdx(order, 2)> dx2() const {
           return this->template dx<2>();
       };

       /*!
        * Calculates the third derivative.
        */
       myspline<T, orderdx(order, 3)> dx3() const {
           return this->template dx<3>();
       };

}; // class myspline
//####################################################### End of defintion of myspline class ############################################################################################

/*!
 * Commutation of spline scalar multiplication operator.
 *
 * @param d Scalar to be multiplied.
 * @param b Spline to be multiplied.
 * @tparam T Datatype of spline and scalar.
 * @tparam order Order of the spline.
 */
template<typename T, size_t order>
inline myspline<T, order> operator*(const T& d, const myspline<T, order> &b) {
    return b * d;
};


/*!
 * The methods in the internal namespace are not supposed to be called from outside this header file.
 */
namespace internal {

/*!
 * Find the overlapping intervals of splines m1 and m2. startindex1 is set to the index of the start of the first overlapping interval in
 * _intervals of m1 and startindex2 to the start of the first overlapping interval in _intervals of m2. nintervals is set to the number
 * of intervals on which m1 and m2 overlap.
 *
 * If there is no overlap startindex1, startindex2 and nintervals are all set to 0.
 *
 * @param m1 First spline.
 * @param m2 Second spline.
 * @param startindex1 Return reference for the index of the first overlapping interval in the intervals vector of spline m1.
 * @param startindex2 Return reference for the index of the first overlapping interval in the intervals vector of spline m2.
 * @param nintervals Return reference for the number of overlapping intervals.
 * @tparam T Datatype of both splines.
 * @tparam order1 Order of spline m1.
 * @tparam order2 Order of spline m2.
 */
template<typename T, size_t order1, size_t order2>
void findOverlappingIntervals(const myspline<T, order1> &m1, const myspline<T, order2> &m2, size_t &startindex1, size_t &startindex2, size_t &nintervals) {
    if (!m1.checkOverlap(m2)) {
        startindex1 = 0;
        startindex2 = 0;
        nintervals = 0;
        return;
    }

    if(m2.getIntervals().front() <= m1.getIntervals().front()) {
        startindex1 = 0;
        startindex2 = std::distance(m2.getIntervals().begin(), std::find(m2.getIntervals().begin(), m2.getIntervals().end(), m1.getIntervals().front()));
        assert(m2.getIntervals().size() >= startindex2 + 2);
    } else {
        startindex2 = 0;
        startindex1 = std::distance(m1.getIntervals().begin(), std::find(m1.getIntervals().begin(), m1.getIntervals().end(), m2.getIntervals().front()));
        assert(m1.getIntervals().size() >= startindex1 + 2);
    }
    nintervals = std::min(m2.getIntervals().size() - startindex2-1, m1.getIntervals().size() - startindex1-1);
};


}; // end namespace myspline::internal

/*!
 * Generates a Bspline of order k-1 at knot i relative to the grid given by knots.
 *
 * @param knots Grido on which to generate the BSplines.
 * @param i Index of the knot at which to generate the BSpline.
 * @tparam T Datatype of the spline.
 * @tparam k Number of the coefficients per interval for the spline (i.e. order of the spline plus one).
 */
template<typename T, size_t k>
myspline<T, k-1> generateBspline(const std::vector<T> &knots, size_t i){
    static_assert(k >= 1, "k has to be at least 1.");
    if constexpr (k == 1) {
        const T& xi = knots.at(i);
        const T& xip1 = knots.at(i+1);
        assert(xip1 > xi);
        std::vector<T> intervals{xi, xip1};
        std::vector<std::array<T, 1>> coefficients{{static_cast<T>(1)}};
        return myspline<T, 0>(std::move(intervals), std::move(coefficients));
    } else {
        myspline<T, k-1> ret;

        const T& xi = knots.at(i);
        const T& xipkm1 = knots.at(i+k-1);
         if (xipkm1 > xi) {
             myspline<T,k-2> spline1 = generateBspline<T, k-1>(knots, i);
             T prefac = static_cast<T>(1)/(xipkm1 -xi);
             spline1 *= prefac; 
             ret += spline1.timesx() - xi * spline1;
        }

        const T& xip1 = knots.at(i+1);
        const T& xipk = knots.at(i+k);
        if (xipk > xip1) {
            myspline<T,k-2> spline2 = generateBspline<T, k-1>(knots, i + 1);
            T prefac = static_cast<T>(1)/(xipk -xip1);
            spline2 *=  prefac;
            ret += xipk * spline2 - spline2.timesx();
        }
        return ret;
    }
};

};
#endif // MYSPLINE_H
