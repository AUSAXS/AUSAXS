#pragma once

#include <Tools.h>
#include <string>

using std::string;

/**
 * @brief \class Reader. 
 *               Virtual super-class for all data file readers. 
 */
class Reader {
  public:
    /**
     * @brief Read the data stored in a file. 
     */
    virtual void read(const string&) {
        print_err("FATAL ERROR: This code should be unreachable.");
        exit(1);
    }
};