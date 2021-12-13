#pragma once

#include <Tools.h>
#include <string>

using std::string;

class Writer {
    public:
        /**
         * @brief Read the file backing this File object. 
         */
        virtual void write(const string&) {
            print_err("FATAL ERROR: This code should be unreachable.");
            exit(1);
        }
};