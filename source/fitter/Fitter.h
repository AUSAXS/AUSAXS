#pragma once

#include <string>
#include <vector>

class Fitter {
public:
    struct Fit {
        std::vector<double> params;
    };

    /**
     * @brief Prepare a new fit for the input file. 
     * @param input the data file to be fitted. 
     */
    Fitter(std::string input) : file(input) {read();}
    virtual ~Fitter() {}
    virtual Fit fit() const = 0;

protected:
    const std::string filename;
    virtual void read() = 0;
};