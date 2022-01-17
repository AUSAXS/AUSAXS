#pragma once

#include <Tools.h>
#include <string>

using std::string;

/**
 * @brief \class Writer. 
 *               Virtual super-class for all data file writers. 
 */
class Writer {
  public:
    /**
     * @brief Write the contents of the backing File to a given path. 
     */
    virtual void write(const string&) {
        print_err("FATAL ERROR: This code should be unreachable.");
        exit(1);
    }
};