/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#ifndef BSPLINE_SUPPORT_SUPPORT_H
#define BSPLINE_SUPPORT_SUPPORT_H

#include <bspline/exceptions/BSplineException.h>
#include <bspline/support/Grid.h>

#include <optional>

/*!
 * The namspace containing the internal representation of the Spline's support
 * and the global grid.
 *
 * @brief Namespace for the spline's grid and support.
 */
namespace bspline::support {
using namespace bspline::exceptions;

/*!
 * Represents the support of a spline as a set of intervals, represented by the
 * the corresponding grid points. It is essentially a view onto the global grid.
 *
 * @brief Represents the spline's support.
 * @tparam T Datatype of the grid and spline.
 */
template <typename T>
class Support final {
 public:
  /*!
   * Represents an interval relative to the global grid represented by _grid.
   * @brief Represents an index relative to the global grid.
   */
  using AbsoluteIndex = size_t;

  /*!
   * Represents an interval relative to the intervals contained in this Support.
   * @brief Represents an index relative to the Support.
   */
  using RelativeIndex = size_t;

  /*!
   * @brief The iterator type.
   */
  using const_iterator = typename Grid<T>::const_iterator;

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

  /*!
   * Checks the validity of this Support and throws if the Support is not in a
   * valid state.
   *
   * @brief Checks whether this Support fulfills the classes invariants.
   * @throws BSplineException if this object is not in a valid state.
   */
  void checkValidity() const {
    const bool isEmpty = _startIndex == 0 && _endIndex == 0;
    const bool containsElement = _endIndex > _startIndex;
    const bool withinBounds = _endIndex <= _grid.size();

    const bool valid = isEmpty || (containsElement && withinBounds);

    if (!valid) {
      throw BSplineException(ErrorCode::INCONSISTENT_DATA);
    }
  }

 public:
  /*!
   * @brief Constructs a Support relative to the global grid grid.
   * @param grid The global grid.
   * @param startIndex The index of the first grid point which is part of the
   * Support.
   * @param endIndex The index of the element behind the last grid point which
   * is part of the Support.
   */
  Support(const Grid<T> &grid, AbsoluteIndex startIndex, AbsoluteIndex endIndex)
      : _grid(grid), _startIndex(startIndex), _endIndex(endIndex) {
    checkValidity();
  };

  /*!
   * @brief Default copy constructor.
   * @param s The Support to copy.
   */
  Support(const Support &s) noexcept = default;

  /*!
   * @brief Default copy assignment operator.
   * @param s The Support to copy.
   */
  Support &operator=(const Support &s) noexcept = default;

  /*!
   * @brief Default destructor.
   */
  ~Support() = default;

  /*!
   * Move constructor. Leaves the moved-from object as an empty Support relative
   * to the same grid.
   *
   * @brief Move constructor.
   * @param s The Support to move.
   */
  Support(Support &&s) noexcept
      : _grid{s._grid}, _startIndex{s._startIndex}, _endIndex{s._endIndex} {
    s._startIndex = 0;
    s._endIndex = 0;
  }

  /*!
   * Move assignment operator. Leaves the moved-from object as an empty Support
   * relative to the same grid.
   *
   * @brief Move assignment operator.
   * @param s The Support to move.
   */
  Support &operator=(Support &&s) noexcept {
    _grid = s._grid;
    _startIndex = s._startIndex;
    _endIndex = s._endIndex;

    s._startIndex = 0;
    s._endIndex = 0;

    return *this;
  }

  /*!
   * @brief Constructs an empty Support relative to the global grid grid.
   * @param grid The global grid.
   */
  static Support<T> createEmpty(const Grid<T> &grid) {
    return Support<T>{grid, 0, 0};
  };

  /*!
   * @brief Constructs a Support representing the complete global grid grid.
   * @param grid The global grid.
   */
  static Support<T> createWholeGrid(const Grid<T> &grid) {
    const auto gridSize = grid.size();
    return Support<T>{grid, 0, gridSize};
  };

  /*!
   * @brief Returns the number of grid points contained in this Support.
   * @returns The number of grid points contained in this Support.
   */
  size_t size() const {
    DURING_TEST_CHECK_VALIDITY();
    return _endIndex - _startIndex;
  };

  /*!
   * Checks whether the number of grid points contained in the Support is zero.
   *
   * @brief Returns whether this Support is empty.
   * @returns True if this Support contains no grid points, false otherwise.
   */
  bool empty() const {
    DURING_TEST_CHECK_VALIDITY();
    return (_startIndex == _endIndex);
  };

  /*!
   * Checks whether the Support contains any intervals. The number of intervals
   * is the number of grid points minus one (size() - 1).
   *
   * @brief Returns whether this Support is empty or point-like.
   * @returns false if the Support is empty or point-like, true otherwise.
   */
  bool containsIntervals() const {
    DURING_TEST_CHECK_VALIDITY();
    return (size() > 1);
  };

  /*!
   *
   * Returns the index relative to this Support from an index relative to the
   * global grid. If the global index does not correspond to a grid point
   * contained in this Support, std::nullopt is returned.
   *
   * @brief Converts an AbsoluteIndex into a RelativeIndex.
   * @param index The AbsoluteIndex referring to an interval on the global grid.
   * @returns The relative index if the corresponding grid point is part of this
   * Support, std::nullopt else.
   */
  std::optional<RelativeIndex> relativeFromAbsolute(AbsoluteIndex index) const {
    DURING_TEST_CHECK_VALIDITY();
    if (index >= _startIndex && index < _endIndex) {
      return index - _startIndex;
    } else {
      return std::nullopt;
    }
  };

  /*!
   * Returns the index relative to this Support from an index relative to the
   * global grid. If the global index does not correspond to an interval
   * contained in this Support, std::nullopt is returned.
   *
   * @brief Converts an AbsoluteIndex into a RelativeIndex.
   * @param index The AbsoluteIndex referring to an interval on the global grid.
   * @returns The relative index if the corresponding interval is part of this
   * Support, std::nullopt else.
   */
  std::optional<RelativeIndex> intervalIndexFromAbsolute(
      AbsoluteIndex index) const {
    DURING_TEST_CHECK_VALIDITY();
    if (index >= _startIndex && index + 1 < _endIndex) {
      return index - _startIndex;
    } else {
      return std::nullopt;
    }
  };

  /*!
   * Returns the index relative to the global grid from an inde relative to this
   * Support.
   *
   * @brief Converts a RelativeIndex into an AbsoluteIndex.
   * @param index The AbsoluteIndex referring to an interval on the global grid.
   * @throws BSplineException If the relative index is out of bounds for this
   * Support.
   * @returns The absolute index corresponding to the relative index.
   */
  AbsoluteIndex absoluteFromRelative(RelativeIndex index) const {
    DURING_TEST_CHECK_VALIDITY();
    if (index >= size()) {
      throw BSplineException(ErrorCode::UNDETERMINED);
    }
    return index + _startIndex;
  };

  /*!
   * @brief Returns the number of intervals represented by this Support.
   * @returns The number of intervals contained in this Support.
   */
  size_t numberOfIntervals() const {
    DURING_TEST_CHECK_VALIDITY();
    const size_t si = size();
    if (si == 0) {
      return 0;
    } else {
      return si - 1;
    }
  };

  /*!
   * @brief Returns the global grid.
   * @returns A reference to the global grid.
   */
  const Grid<T> &getGrid() const {
    DURING_TEST_CHECK_VALIDITY();
    return _grid;
  };

  /*!
   * @brief Returns the _startIndex.
   * @returns The start index.
   */
  AbsoluteIndex getStartIndex() const {
    DURING_TEST_CHECK_VALIDITY();
    return _startIndex;
  };

  /*!
   * @brief Returns the _endIndex.
   * @returns The end index.
   */
  AbsoluteIndex getEndIndex() const {
    DURING_TEST_CHECK_VALIDITY();
    return _endIndex;
  }

  /*!
   * Allows access to the grid points contained in the Support. Performs no
   * bounds checks.
   *
   * @brief Returns the indexed element.
   * @param index Index of the element.
   * @returns A reference to the element represented by index. Bounds are not
   * checked.
   */
  const T &operator[](RelativeIndex index) const {
    DURING_TEST_CHECK_VALIDITY();
    return _grid[_startIndex + index];
  };

  /*!
   * Allows access to the grid points contained in the Support. Checks bounds
   * and throws exception in case of out-of-bounds access.
   *
   * @brief Returns the indexed element.
   * @param index Index of the element.
   * @throws BSplineException If the access is out of bounds for this Support.
   * @returns A reference to the element represented by index. Bounds are not
   * checked.
   */
  const T &at(RelativeIndex index) const {
    DURING_TEST_CHECK_VALIDITY();
    if (_startIndex + index >= _endIndex) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _grid.at(_startIndex + index);
  };

  /*!
   * Returns a reference to the first grid point that is part of the Support.
   *
   * @brief Returns the begin of this Support.
   * @throws BSplineException If this Support is empty.
   * @returns A reference to the first grid point contained in this Support.
   */
  const T &front() const {
    DURING_TEST_CHECK_VALIDITY();
    if (empty()) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _grid[_startIndex];
  };

  /*!
   * Returns a reference to the last grid point that is part of the Support.
   *
   * @brief Returns the end of this Support.
   * @throws BSplineException If this Support is empty.
   * @returns A reference to the last grid point contained in this Support.

   */
  const T &back() const {
    DURING_TEST_CHECK_VALIDITY();
    if (empty()) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _grid[_endIndex - 1];
  };

  /*!
   * @brief Returns the begin iterator of the Support.
   * @returns An iterator to the first element.
   */
  const_iterator begin() const {
    DURING_TEST_CHECK_VALIDITY();
    return _grid.begin() + _startIndex;
  };

  /*!
   * @brief Returns the end iterator of the Support.
   * @returns An iterator pointing behind the last element.
   */
  const_iterator end() const {
    DURING_TEST_CHECK_VALIDITY();
    return _grid.begin() + _endIndex;
  };

  /*!
   * Checks whether the global Grids, the two Supports are defined on, are
   * logically equivalent.
   *
   * @brief Checks whether two Supports are defined on the same Grids.
   * @param s Support to check against.
   * @returns True if the Support s is defined on a logically equivalent grid,
   * false otherwise.
   */
  bool hasSameGrid(const Support &s) const {
    DURING_TEST_CHECK_VALIDITY();
    return _grid == s._grid;
  };

  /*!
   * Compares two Supports for equality. For two Supports two be equal, they
   * have to be defined on the same grid, must represent the same subset of the
   * number line.
   *
   * @brief Logically compares to Supports.
   * @param s Support to compare against.
   * @returns True it the two Supports are logically equivalent, false
   * otherwise.
   */
  bool operator==(const Support &s) const {
    DURING_TEST_CHECK_VALIDITY();
    DURING_TEST_CHECK_VALIDITY_OF(s);
    return hasSameGrid(s) &&
           ((_startIndex == s._startIndex && _endIndex == s._endIndex) ||
            (empty() && s.empty()));
  };

  /*!
   * Compares two Supports for inequality. For two Supports two be equal, they
   * have to be defined on the same grid, must represent the same subset of the
   * number line.
   *
   * @brief Logically compares to Supports.
   * @param s Support to compare against.
   * @returns False it the two Supports are logically equivalent, true
   * otherwise.
   */
  bool operator!=(const Support &s) const {
    DURING_TEST_CHECK_VALIDITY();
    return !(*this == s);
  };

  /*!
   * Calculates the union of this Support with the Support s. This is not
   * strictly the set-theoretical union (if the two Supports do not overlap),
   * but a Support representing one contiguous bit of the number line containing
   * both Supports.
   *
   * @brief Calculates the union of the two Supports.
   * @param s The Support to calculate the union with.
   * @throws BSplineException If the two Supports are defined on (logically)
   * different grids.
   * @returns The Support representing the union of the two Supports.
   */
  Support calcUnion(const Support &s) const {
    DURING_TEST_CHECK_VALIDITY();
    DURING_TEST_CHECK_VALIDITY_OF(s);
    if (!hasSameGrid(s)) {
      throw BSplineException(ErrorCode::DIFFERING_GRIDS);
    }

    const bool thisEmpty = empty();
    const bool sEmpty = s.empty();

    if (thisEmpty && sEmpty) {
      // Both Supports are empty, return empty Support.
      return createEmpty(_grid);
    } else if (thisEmpty && !sEmpty) {
      return s;
    } else if (!thisEmpty && sEmpty) {
      return *this;
    }

    const size_t newStartIndex = std::min(_startIndex, s._startIndex);
    const size_t newEndIndex = std::max(_endIndex, s._endIndex);
    return Support(_grid, newStartIndex, newEndIndex);
  };

  /*!
   * Calculates the intersection of the two Supports.
   *
   * @brief Calculates the intersection of the two Supports.
   * @param s The Support to calculate the intersection with.
   * @throws BSplineException If the two Supports are defined on (logically)
   * different grids.
   * @returns The Support representing the intersection of the two Supports.
   */
  Support calcIntersection(const Support &s) const {
    DURING_TEST_CHECK_VALIDITY();
    DURING_TEST_CHECK_VALIDITY_OF(s);
    if (!hasSameGrid(s)) {
      throw BSplineException(ErrorCode::DIFFERING_GRIDS);
    }

    const size_t newStartIndex = std::max(_startIndex, s._startIndex);
    const size_t newEndIndex = std::min(_endIndex, s._endIndex);

    if (newStartIndex >= newEndIndex) {
      // No overlap, return empty Support.
      return createEmpty(_grid);
    } else {
      return Support(_grid, newStartIndex, newEndIndex);
    }
  };
};  // end class Support
}  // namespace bspline::support
#endif  // BSPLINE_SUPPORT_SUPPORT_H
