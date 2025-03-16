/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#ifndef BSPLINE_CONCEPTS_H
#define BSPLINE_CONCEPTS_H

#include <concepts>
#include <type_traits>

/*!
 * @brief Main namespace for this library.
 */
namespace bspline {
#ifndef BSPLINE_DOXYGEN_IGNORE
namespace internal {

/*!
 * Template class to tell the library to treat a given type as a real value.
 *
 * @tparam T The type to be defined as real.
 */
template <typename T>
class is_registered_real : public std::false_type {};
}  // namespace internal
#endif

/*!
 * Concept that makes the library treat a type as a real value.
 */
template <typename T>
concept Real = std::floating_point<T> || internal::is_registered_real<T>::value;

/*!
 * Concept that makes the library treat a type as a scalar.
 */
template <typename T>
concept Scalar = Real<T> || std::integral<T>;

}  // namespace bspline

#define REGISTER_BSPLINE_REAL(Type)                          \
  namespace bspline::internal {                              \
  template <>                                                \
  class is_registered_real<Type> : public std::true_type {}; \
  }                                                          \
  static_assert(true)

#endif  // BSPLINE_CONCEPTS_H
