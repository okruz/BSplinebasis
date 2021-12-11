#ifndef OKRUZ_BSPLINE_SUPPORT_GRID_H
#define OKRUZ_BSPLINE_SUPPORT_GRID_H
#include <okruz/bspline/exceptions/BSplineException.h>
#include <okruz/bspline/internal/misc.h>

#include <algorithm>
#include <memory>
#include <vector>

/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

namespace okruz::bspline::support {
using namespace okruz::bspline::exceptions;

/*!
 * Represents a global grid.
 *
 * @tparam T The datatype of the grid elements.
 */
template <typename T>
class Grid {
 private:
  /*! The gridpoints. */
  std::shared_ptr<const std::vector<T>> _data;

  bool isSteadilyIncreasing() const {
    for (size_t i = 1; i < _data->size(); i++) {
      if ((*_data)[i - 1] >= (*_data)[i]) {
        return false;
      }
    }
    return true;
  }

 public:
  /*!
   * Iterator type.
   */
  using const_iterator = typename std::vector<T>::const_iterator;

  Grid() = delete;

  /*!
   * Constructs a grid from a set of begin and end iterators.
   *
   * @param begin The iterator referencing the first element to be copied into
   * the grid.
   * @param end The iterator referencing the element behind the last element to
   * be copied into the grid.
   * @tparam Iter The type of the two iterators.
   */
  template <typename Iter>
  Grid(Iter begin, Iter end)
      : Grid(std::make_shared<const std::vector<T>>(begin, end)) {}

  /*!
   * Constructs a grid from a std::vector.
   *
   * @param v The input vector.
   */
  explicit Grid(std::vector<T> v)
      : Grid(std::make_shared<const std::vector<T>>(std::move(v))){};

  /*!
   * Constructs a grid from a std::initializer_list.
   *
   * @param v The input initializer_list.
   */
  explicit Grid(const std::initializer_list<T> &v) : Grid(v.begin(), v.end()){};

  /*!
   * Constructs a grid from a std::shared_ptr<const std::vector<T>>.
   *
   * @param data A shared pointer to the grid elements.
   */
  explicit Grid(std::shared_ptr<const std::vector<T>> data)
      : _data(std::move(data)) {
    if (!_data || !isSteadilyIncreasing()) {
      throw BSplineException(ErrorCode::INCONSISTENT_DATA,
                             "The grid points are not steadily increasing.");
    }
  };

  /*!
   * Comparison operator.
   *
   * @param g The grid to compare this grid with.
   * @returns Returns true if the grids represent the same logical grid.
   */
  bool operator==(const Grid &g) const {
    if (_data == g._data) {
      // Fast path: Compare pointers.
      return true;
    }

    // Slow path: Compare logically.
    return (*_data) == (*g._data);
  }

  /*!
   * Returns the number of elements of the grid.
   */
  size_t size() const { return _data->size(); };

  /*!
   * Returns a shared pointer to the elements of this grid.
   */
  std::shared_ptr<const std::vector<T>> getData() const { return _data; };

  /*!
   * Checks whether this spline holds no elements.
   *
   * @returns Returns true if this grid holds no element.
   */
  bool empty() const { return _data->empty(); };

  /*!
   * Returns a reference to the ith element of the grid. Performs no bounds
   * checks.
   *
   * @param i The index of the element to be returned.
   * @returns A reference to the ith element.
   */
  const T &operator[](size_t i) const { return (*_data)[i]; };

  /*!
   * Returns a reference to the ith element of the grid. Checks the bounds.
   *
   * @param i The index of the element to be returned.
   * @returns A reference to the ith element.
   */
  const T &at(size_t i) const {
    if (i > size()) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return (*_data)[i];
  };

  /*!
   * Returns a reference to the first element of the grid.
   *
   * @returns A reference to the first element.
   */
  const T &front() const {
    if (empty()) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _data->front();
  };

  /*!
   * Returns a reference to the last element of the grid.
   *
   * @returns A reference to the last element.
   */
  const T &back() const {
    if (empty()) {
      throw BSplineException(ErrorCode::INVALID_ACCESS);
    }
    return _data->back();
  };

  /*!
   * Returns the begin iterator of the grid.
   *
   * @returns An iterator to the first element.
   */
  const_iterator begin() const { return _data->begin(); };

  /*!
   * Returns the end iterator of the grid..
   *
   * @returns An iterator pointing behind the last element.
   */
  const_iterator end() const { return _data->end(); };

  /*!
   * Returns the index corresponding to the element x.
   *
   * @param x The element to be searched for.
   */
  size_t findElement(const T &x) const {
    const auto beginIt = begin();
    const auto endIt = end();
    const auto elementIt = std::lower_bound(beginIt, endIt, x);

    if (elementIt == endIt || *elementIt != x) {
      // Element was not found
      throw BSplineException(ErrorCode::INCONSISTENT_DATA);
    } else {
      return std::distance(beginIt, elementIt);
    }
  };
};
}  // namespace okruz::bspline::support
#endif  // OKRUZ_BSPLINE_SUPPORT_GRID_H
