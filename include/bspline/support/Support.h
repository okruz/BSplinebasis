#ifndef BSPLINE_SUPPORT_SUPPORT_H
#define BSPLINE_SUPPORT_SUPPORT_H

#include <bspline/exceptions/BSplineException.h>
#include <bspline/support/Grid.h>

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
namespace bspline::support {
using namespace bspline::exceptions;

/*!
 * Represents the support of a spline as a set of intervals, represented by the
 * the corresponding grid points. It is essentially a view onto the global grid.
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
   */
  static Support<T> createEmpty(Grid<T> grid) {
    return Support<T>{std::move(grid), 0, 0};
  };

  /*!
   * Constructs a support representing the complete global grid grid.
   *
   * @param grid The global grid.
   */
  static Support<T> createWholeGrid(Grid<T> grid) {
    const auto gridSize = grid.size();
    return Support<T>{std::move(grid), 0, gridSize};
  };

  /*!
   * Returns the number of grid points contained in this support.
   *
   * @returns The number of grid points contained in this support.
   */
  size_t size() const { return _endIndex - _startIndex; };

  /*!
   * Checks whether the number of grid points contained in the support is zero.
   *
   * @returns True if this support contains no grid points, false otherwise.
   */
  bool empty() const { return (_startIndex == _endIndex); };

  /*!
   * Checks whether the support contains any intervals. The number of intervals
   * is the number of grid points minus one (size() - 1).
   *
   * @returns false if the support is empty or point-like, true otherwise.
   */
  bool containsIntervals() const { return (size() > 1); };

  /*!
   *
   * Returns the index relative to this support from an index relative to the
   * global grid. If the global index does not correspond to a grid point
   * contained in this support, std::nullopt is returned.
   *
   * @param index The AbsoluteIndex referring to an interval on the global grid.
   * @returns The relative index if the corresponding grid point is part of this
   * support, std::nullopt else.
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
   * @returns The relative index if the corresponding interval is part of this
   * support, std::nullopt else.
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
   * @throws BSplineException If the relative index is out of bounds for this
   * support.
   * @returns The absolute index corresponding to the relative index.
   */
  AbsoluteIndex absoluteFromRelative(RelativeIndex index) const {
    if (index >= size()) {
      throw BSplineException(ErrorCode::UNDETERMINED);
    }
    return index + _startIndex;
  };

  /*!
   * Returns the number of intervals represented by this support.
   *
   * @returns The number of intervals contained in this support.
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
   *
   * @returns A reference to the global grid.
   */
  const Grid<T> &getGrid() const { return _grid; };

  /*!
   * Returns the _startIndex.
   *
   * @returns The start index.
   */
  AbsoluteIndex getStartIndex() const { return _startIndex; };

  /*!
   * Returns the _endIndex.
   *
   * @returns The end index.
   */
  AbsoluteIndex getEndIndex() const { return _endIndex; }

  /*!
   * Allows access to the grid points contained in the support. Performs no
   * bounds checks.
   *
   * @param index Index of the element.
   * @returns A reference to the element represented by index. Bounds are not
   * checked.
   */
  const T &operator[](RelativeIndex index) const {
    return _grid[_startIndex + index];
  };

  /*!
   * Allows access to the grid points contained in the support. Checks bounds
   * and throws exception in case of out-of-bounds access.
   *
   * @param index Index of the element.
   * @throws BSplineException If the access is out of bounds for this support.
   * @returns A reference to the element represented by index. Bounds are not
   * checked.
   */
  const T &at(RelativeIndex index) const {
    if (_startIndex + index >= _endIndex) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _grid.at(_startIndex + index);
  };

  /*!
   * Returns a reference to the first grid point that is part of the support.
   *
   * @throws BSplineException If this support is empty.
   * @returns A reference to the first grid point contained in this support.
   */
  const T &front() const {
    if (empty()) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _grid[_startIndex];
  };

  /*!
   * Returns a reference to the last grid point that is part of the support.
   *
   * @throws BSplineException If this support is empty.
   * @returns A reference to the last grid point contained in this support.

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
   * @returns True if the support s is defined on a logically equivalent grid,
   * false otherwise.
   */
  bool hasSameGrid(const Support &s) const { return _grid == s._grid; };

  /*!
   * Compares two supports for equality. For two supports two be equal, they
   * have to be defined on the same grid, must represent the same subset of the
   * number line.
   *
   * @param s Support to compare against.
   * @returns True it the two supports are logically equivalent, false
   * otherwise.
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
   * @throws BSplineException If the two supports are defined on (logically)
   * different grids.
   * @returns The support representing the union of the two supports.
   */
  Support calcUnion(const Support &s) const {
    if (!hasSameGrid(s)) {
      throw BSplineException(ErrorCode::DIFFERING_GRIDS);
    }
    const bool thisEmpty = empty();
    const bool sEmpty = s.empty();
    if (thisEmpty && sEmpty)
      return createEmpty(
          _grid);  // Both Supports are empty, return empty Support
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
   * @throws BSplineException If the two supports are defined on (logically)
   * different grids.
   * @returns The support representing the intersection of the two supports.
   */
  Support calcIntersection(const Support &s) const {
    if (!hasSameGrid(s)) {
      throw BSplineException(ErrorCode::DIFFERING_GRIDS);
    }
    const size_t newStartIndex = std::max(_startIndex, s._startIndex);
    const size_t newEndIndex = std::min(_endIndex, s._endIndex);
    if (newStartIndex >= newEndIndex)
      return createEmpty(_grid);  // no overlap, return empty Support
    else
      return Support(_grid, newStartIndex, newEndIndex);
  };
};  // end class Support
}  // namespace bspline::support
#endif  // BSPLINE_SUPPORT_SUPPORT_H
