#include "Fitter.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

using std::string, std::vector;

class IntensityFitter : public Fitter {
public: 
    using Fitter::Fitter;
    ~IntensityFitter() override {}
    Fit fit() const override {}

private: 
    void read() override {
        // check if file was succesfully opened
        std::ifstream input(filename);
        if (!input.is_open()) {throw std::ios_base::failure("Error in IntensityFitter::read: Could not open file \"" + filename + "\"");}
        
        string line; // placeholder for the current line
        while(getline(input, line)) {
            vector<string> tokens;
            boost::split(tokens, line, boost::is_any_of(" ,\t")); // spaces, commas, and tabs can all be used as separators

            // determine if we are in some sort of header
            if (tokens.size() < 3 || tokens.size() > 4) {continue;} // too many separators
            for (const auto& e : tokens) { // check if they are numbers
                if (!e.empty() && e.find_first_not_of("0123456789,") != string::npos) {continue;}
            }

            // now we are most likely beyond any headers
            double q, I, sigma;
            q = std::stod(tokens[0]); // we know for sure that the strings are convertible to numbers (boost check)
            I = std::stod(tokens[1]);
            sigma = std::stod(tokens[2]);

            std::cout << q << ", " << I << ", " << sigma << std::endl; 
        }
    }
};