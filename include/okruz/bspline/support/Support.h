#ifndef OKRUZ_BSPLINE_SUPPORT_SUPPORT_H
#define OKRUZ_BSPLINE_SUPPORT_SUPPORT_H

#include <okruz/bspline/exceptions/BSplineException.h>
#include <okruz/bspline/support/Grid.h>

#include <optional>

/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

/*!
 * The namspace containing the internal representation of the Spline's support
 * and the global grid.
 */
namespace okruz::bspline::support {
using namespace okruz::bspline::exceptions;

/*!
 * This enum allows to signal to some constructors whether
 * an empty support should be constructed or one representing the
 * whole grid.
 */
enum class Construction { EMPTY, WHOLE_GRID };

/*!
 * Represents the support of a spline as a number of gridpoints.
 * It is hence essentially a view onto the global grid. The pointer to the
 * global grid allows for a speedy comparison whether splines are defined on the
 * same global grid and methods that assume they are can be used.
 *
 * @tparam T Datatype of the grid and spline.
 */
template <typename T>
class Support {
 public:  // Type definitions.
  /*!
   * Represents an interval relative to the global grid represented by _grid.
   */
  using AbsoluteIndex = size_t;

  /*!
   * Represents an interval relative to the intervals contained in this support.
   */
  using RelativeIndex = size_t;

  /*!
   * Iterator type.
   */
  using const_iterator = typename std::vector<T>::const_iterator;

 private:
  /*! Represents the global grid. */
  Grid<T> _grid;
  /*! Represents the begin of the Support. */
  AbsoluteIndex _startIndex;
  /*!
   * Represents the end of the Support. Points to the element behind the last
   * element of the Support.
   */
  AbsoluteIndex _endIndex;

 public:
  /*!
   * Constructs a support relative to the global grid grid.
   *
   * @param grid The global grid.
   * @param startIndex The index of the first grid point which is part of the
   * support.
   * @param endIndex The index of the element behind the last grid point which
   * is part of the support.
   */
  Support(Grid<T> grid, AbsoluteIndex startIndex, AbsoluteIndex endIndex)
      : _grid(std::move(grid)), _startIndex(startIndex), _endIndex(endIndex) {
    if (_startIndex > _endIndex || _endIndex > _grid.size()) {
      throw BSplineException(ErrorCode::INCONSISTENT_DATA);
    }
  };

  /*!
   * Constructs an empty support relative to the global grid grid.
   *
   * @param grid The global grid.
   * @param constr The construction. If this parameter is not provided, an empty
   * support will be constructed, if it is set to Construction::WHOLE_GRID a
   * support representing the extend of the whole grid will be constructed.
   */
  Support(Grid<T> grid, Construction constr = Construction::EMPTY)
      : _grid(std::move(grid)),
        _startIndex(0),
        _endIndex((constr == Construction::EMPTY) ? 0 : _grid.size()){};

  /*!
   * Returns the number of grid points contained in the support.
   */
  size_t size() const { return _endIndex - _startIndex; };

  /*!
   * Checks whether the number of grid points contained in the support is zero.
   */
  bool empty() const { return (_startIndex == _endIndex); };

  /*!
   * Checks whether the support contains any intervals. Returns true if the
   * support is empty or point-like. The number of intervals is the number of
   * grid points minus one (size() - 1).
   */
  bool containsIntervals() const { return (size() > 1); };

  /*!
   *
   * Returns the index relative to this support from an index relative to the
   * global grid. If the global index does not correspond to a grid point
   * contained in this support, std::nullopt is returned.
   *
   * @param index The AbsoluteIndex referring to an interval on the global grid.
   */
  std::optional<RelativeIndex> relativeFromAbsolute(AbsoluteIndex index) const {
    if (index >= _startIndex && index < _endIndex)
      return index - _startIndex;
    else
      return std::nullopt;
  };

  /*!
   * Returns the index relative to this support from an index relative to the
   * global grid. If the global index does not correspond to an interval
   * contained in this support, std::nullopt is returned.
   *
   * @param index The AbsoluteIndex referring to an interval on the global grid.
   */
  std::optional<RelativeIndex> intervalIndexFromAbsolute(
      AbsoluteIndex index) const {
    if (index >= _startIndex && index + 1 < _endIndex)
      return index - _startIndex;
    else
      return std::nullopt;
  };

  /*!
   * Returns the index relative to the global grid from an inde relative to this
   * support.
   *
   * @param index The AbsoluteIndex referring to an interval on the global grid.
   */
  AbsoluteIndex absoluteFromRelative(RelativeIndex index) const {
    if (index >= size()) {
      throw BSplineException(ErrorCode::UNDETERMINED);
    }
    return index + _startIndex;
  };

  /*!
   * Returns the number of intervals represented by this support.
   */
  size_t numberOfIntervals() const {
    const size_t si = size();
    if (si == 0)
      return 0;
    else
      return si - 1;
  };

  /*!
   * Returns the global grid.
   */
  Grid<T> getGrid() const { return _grid; };

  /*!
   * Returns the _startIndex.
   */
  AbsoluteIndex getStartIndex() const { return _startIndex; };

  /*!
   * Returns the _endIndex.
   */
  AbsoluteIndex getEndIndex() const { return _endIndex; }

  /*!
   * Allows access to the grid points contained in the support. Performs no
   * bounds checks.
   *
   * @param index Index of the element.
   */
  const T &operator[](RelativeIndex index) const {
    return _grid[_startIndex + index];
  };

  /*!
   * Allows access to the grid points contained in the support. Checks bounds
   * and throws exception in case of out-of-bounds access.
   *
   * @param index Index of the element.
   */
  const T &at(RelativeIndex index) const {
    if (_startIndex + index >= _endIndex) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _grid.at(_startIndex + index);
  };

  /*!
   * Returns a reference to the first grid point that is part of the support.
   */
  const T &front() const {
    if (empty()) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _grid[_startIndex];
  };

  /*!
   * Returns a reference to the last grid point that is part of the support.
   */
  const T &back() const {
    if (empty()) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _grid[_endIndex - 1];
  };

  /*!
   * Returns the begin iterator of the support.
   *
   * @returns An iterator to the first element.
   */
  const_iterator begin() const { return _grid.begin() + _startIndex; };

  /*!
   * Returns the end iterator of the support..
   *
   * @returns An iterator pointing behind the last element.
   */
  const_iterator end() const { return _grid.begin() + _endIndex; };

  /*!
   * Checks whether the global grids, the two supports are defined on, are
   * logically equivalent.
   *
   * @param s Support to check against.
   */
  bool hasSameGrid(const Support &s) const { return _grid == s._grid; };

  /*!
   * Compares two supports for equality. For two supports two be equal, they
   * have to be defined on the same grid, must represent the same subset of the
   * number line.
   *
   * @param s Support to compare against.
   */
  bool operator==(const Support &s) const {
    return hasSameGrid(s) &&
           ((_startIndex == s._startIndex && _endIndex == s._endIndex) ||
            (empty() && s.empty()));
  };

  /*!
   * Calculates the union of this support with the support s. This is not
   * strictly the set-theoretical union (if the two supports do not overlap),
   * but a support representing one contiguous bit of the number line containing
   * both supports.
   *
   * @param s The Support to calculate the union with.
   */
  Support calcUnion(const Support &s) const {
    if (!hasSameGrid(s)) {
      throw BSplineException(ErrorCode::DIFFERING_GRIDS);
    }
    const bool thisEmpty = empty();
    const bool sEmpty = s.empty();
    if (thisEmpty && sEmpty)
      return Support(_grid);  // Both Supports are empty, return empty Support
    else if (thisEmpty && !sEmpty)
      return s;
    else if (!thisEmpty && sEmpty)
      return *this;
    const size_t newStartIndex = std::min(_startIndex, s._startIndex);
    const size_t newEndIndex = std::max(_endIndex, s._endIndex);
    return Support(_grid, newStartIndex, newEndIndex);
  };

  /*!
   * Calculates the intersection of the two Supports.
   *
   * @param s The Support to calculate the intersection with.
   */
  Support calcIntersection(const Support &s) const {
    if (!hasSameGrid(s)) {
      throw BSplineException(ErrorCode::DIFFERING_GRIDS);
    }
    const size_t newStartIndex = std::max(_startIndex, s._startIndex);
    const size_t newEndIndex = std::min(_endIndex, s._endIndex);
    if (newStartIndex >= newEndIndex)
      return Support(_grid);  // no overlap, return empty Support
    else
      return Support(_grid, newStartIndex, newEndIndex);
  };
};  // end class Support
}  // namespace okruz::bspline::support
#endif  // OKRUZ_BSPLINE_SUPPORT_SUPPORT_H
