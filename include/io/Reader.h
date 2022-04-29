#pragma once

#include <string>

/**
 * @brief \class Reader. 
 *               Virtual super-class for all data file readers. 
 */
class Reader {
  public:
    /**
     * @brief Read the data stored in a file. 
     */
    virtual void read(std::string);
};