#pragma once

#include <string>

using std::string;

/**
 * @brief This namespace contains all custom exceptions for this project. 
 */
namespace except {
    struct base : public std::exception {
        base(const char* msg) : msg(msg) {}
        base(const string msg) : msg(msg) {}
        const char* what() const noexcept {return msg.data();}
        const string msg;
    };

    // Invalid call order. A method depends on another before it can be run. Used for fits (a fit must be made before a plot can).
    struct bad_order : public base {using base::base;};

    struct invalid_extension : public base {using base::base;};

    // Invalid argument. Used whenever a check on the arguments is made. 
    struct invalid_argument : public base {using base::base;};

    // An atom is placed out of bounds. Used in the Grid class. 
    struct out_of_bounds : public base {using base::base;};

    // Invalid operation. Used in the Grid class. 
    struct invalid_operation : public base {using base::base;};

    // Unknown string argument. Used in a few different places dealing with user-typed string inputs. 
    struct unknown_argument : public base {using base::base;};

    // Parse error. Used in almost all classes dealing with file inputs with a strict format. 
    struct parse_error : public base {using base::base;};

    // Disabled error. Used when an inherited method is disabled for some reason.  
    struct disabled : public base {using base::base;};

    // Size error. Used when something is wrong with sizes. 
    struct size_error : public base {using base::base;};

    // IO error. Used when something is wrong with reading/writing files. 
    struct io_error : public base {using base::base;};

    // Null-pointer error. Used when a pointer has not been initialized yet. 
    struct nullptr_error : public base {using base::base;};

    // Unexpected error. Used whenever we really did not expect something to go wrong, but it did. 
    struct unexpected : public base {using base::base;};

    // Map error. Used when something is wrong with a map. Used for the constants. 
    struct map_error : public base {using base::base;};
}