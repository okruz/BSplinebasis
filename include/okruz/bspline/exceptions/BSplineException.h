#ifndef OKRUZ_BSPLINE_BSPLINEEXCEPTION_H
#define OKRUZ_BSPLINE_BSPLINEEXCEPTION_H
#include <sstream>
#include <string>

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

/*!
 * The namespace containing the exceptions and errorCodes.
 */
namespace okruz::bspline::exceptions {

/*!
 * The error codes, which may be expected.
 */
enum class ErrorCode {
  DIFFERING_GRIDS,
  INCONSISTENT_DATA,
  INVALID_ACCESS,
  UNDETERMINED
};

/*!
 * Returns the default error message for an error code.
 *
 * @param errorCode The errorCode for which the default message is requested.
 */
std::string getErrorMessage(ErrorCode errorCode) {
  switch (errorCode) {
  case ErrorCode::DIFFERING_GRIDS:
    return "The requested operation is not implemented for splines defined on "
           "different grids.";
  case ErrorCode::INCONSISTENT_DATA:
    return "The data provided is inconsitent.";
  case ErrorCode::INVALID_ACCESS:
    return "Attempted access of nonexistent data.";
  case ErrorCode::UNDETERMINED:
    return "The cause of the error is undetermined.";
  default:
    return "Errorcode unknown.";
  }
};

/*!
 * Returns the errorCode name.
 *
 * @param errorCode The errorCode for which the name is requested.
 */
std::string getErrorCodeName(ErrorCode errorCode) {
  switch (errorCode) {
  case ErrorCode::DIFFERING_GRIDS:
    return "DIFFERING_GRIDS";
  case ErrorCode::INCONSISTENT_DATA:
    return "INCONSISTENT_DATA";
  case ErrorCode::INVALID_ACCESS:
    return "INVALID_ACCESS";
  case ErrorCode::UNDETERMINED:
    return "UNDETERMINED";
  default:
    return "UNKNOWN_ERRORCODE";
  }
};

/*!
 * The main exception class.
 */
class BSplineException : public std::exception {
private:
  /*! The error code. */
  ErrorCode _errorCode;
  /*! A custom defined or defaulted error message. */
  std::string _message;

public:
  /*!
   * Use the default error message for the error code.
   *
   * @param errorCode The errorCode.
   */
  BSplineException(ErrorCode errorCode)
      : _errorCode(errorCode), _message(getErrorMessage(errorCode)){};

  /*!
   * Use the default error message for the error code.
   *
   * @param errorCode The errorCode.
   */

  BSplineException(ErrorCode errorCode, std::string message)
      : _errorCode(errorCode), _message(std::move(message)){};

  /*!
   * Returns a string representaiton of the exception.
   */
  virtual const char *what() const noexcept override {
    std::stringstream ret;
    ret << "BSplineException (code: " << getErrorCodeName(_errorCode)
        << "): " << _message;
    return ret.str().c_str();
  };
};

};     // namespace okruz::bspline::exceptions
#endif // OKRUZ_BSPLINE_BSPLINEEXCEPTION_H
