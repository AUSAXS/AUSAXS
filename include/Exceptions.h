#pragma once

/**
 * @brief This namespace contains all custom exceptions for this project. 
 */
namespace except {
    // Invalid call order. A method depends on another before it can be run. Used for fits (a fit must be made before a plot can).
    struct bad_order : public std::exception {
        bad_order(const char* msg) : msg(msg) {}
        bad_order(const string msg) : msg(msg) {}
        const char* what() const throw() {return msg.data();}
        const string msg;
    };

    // Invalid argument. Used whenever a check on the arguments is made. 
    struct invalid_argument : public std::exception {
        invalid_argument(const char* msg) : msg(msg) {}
        invalid_argument(const string msg) : msg(msg) {}
        const char* what() const throw() {return msg.data();}
        const string msg;
    };

    // An atom is placed out of bounds. Used in the Grid class. 
    struct out_of_bounds : public std::exception {
        out_of_bounds(const char* msg) : msg(msg) {}
        out_of_bounds(const string msg) : msg(msg) {}
        const char* what() const throw() {return msg.data();}
        const string msg;
    };

    // Invalid operation. Used in the Grid class. 
    struct invalid_operation : public std::exception {
        invalid_operation(const char* msg) : msg(msg) {}
        invalid_operation(const string msg) : msg(msg) {}
        const char* what() const throw() {return msg.data();}
        const string msg;
    };

    // Unknown string argument. Used in a few different places dealing with user-typed string inputs. 
    struct unknown_argument : public std::exception {
        unknown_argument(const char* msg) : msg(msg) {}
        unknown_argument(const string msg) : msg(msg) {}
        const char* what() const throw() {return msg.data();}
        const string msg;
    };

    // Parse error. Used in almost all classes dealing with file inputs with a strict format. 
    struct parse_error : public std::exception {
        parse_error(const char* msg) : msg(msg) {}
        parse_error(const string msg) : msg(msg) {}
        const char* what() const throw() {return msg.data();}
        const string msg;
    };
}