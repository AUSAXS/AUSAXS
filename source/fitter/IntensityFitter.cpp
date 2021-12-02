#include "Fitter.h"
#include "math/CubicSpline.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <tuple>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

using std::string, std::vector;

class IntensityFitter : public Fitter {
public: 
    /**
     * @brief Prepare a fit of the measured values in @a input to the model described by @a q and @a I.
     * @param input the path to the file containing the measured values. 
     * @param q the model q values.
     * @param I the model I values. 
     */
    IntensityFitter(string input, vector<double>& q, vector<double>& I) {setup(input, q, I);}
    ~IntensityFitter() override {}

    Fit fit() const override {
        ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
        auto f = std::bind(&IntensityFitter::chi2, this, std::placeholders::_1);
        ROOT::Math::Functor functor(f, 1); // declare the function to be minimized and its number of parameters
        minimizer->SetFunction(functor);
        minimizer->SetVariable(0, "c", 0, 1e-4);
        minimizer->Minimize();
        const double* res = minimizer->X();

        std::map<string, double> pars = {{"c", res[0]}};
        return Fit(pars);
    }

private: 
    vector<double> qo; // observed q values
    vector<double> Io; // observed I values
    vector<double> Im; // model I values
    vector<double> sigma; // error in Io

    /**
     * @brief Calculate chi2 for a given choice of parameters @a params.
     */
    double chi2(const double* params) const {
        double c = params[0];
        double chi = 0;
        for (int i = 0; i < qo.size(); ++i) {
            chi += pow((Io[i] - c*Im[i])/sigma[i], 2);
        }
        return chi;
    }

    /**
     * @brief Prepare this class for fitting.
     * @param file measured values to compare the model against.
     * @param x model q values.
     * @param y model I values.
     */
    void setup(string file, vector<double>& x, vector<double>& y) {
        std::tie(qo, Io, sigma) = read(file); // observed values

        // since our model is not analytically derived, we cannot evaluate it at any random point we'd like
        // thus we have to interpolate 
        Im = vector<double>(qo.size()); // model values
        CubicSpline s(x, y);
        for (int i = 0; i < qo.size(); ++i) {
            Im[i] = s.spline(qo[i]);
        }
    }

    std::tuple<vector<double>, vector<double>, vector<double>> read(string file) const {
        // check if file was succesfully opened
        std::ifstream input(file);
        if (!input.is_open()) {throw std::ios_base::failure("Error in IntensityFitter::read: Could not open file \"" + file + "\"");}

        vector<double> q;
        vector<double> I;
        vector<double> sigma;
        string line; // placeholder for the current line
        while(getline(input, line)) {
            if (line[0] == ' ') {line = line.substr(1);} // fix leading space
            vector<string> tokens;
            boost::split(tokens, line, boost::is_any_of(" ,\t")); // spaces, commas, and tabs can all be used as separators (but not a mix of them)

            // determine if we are in some sort of header
            if (tokens.size() < 3 || tokens.size() > 4) {continue;} // too many separators
            bool skip = false;
            for (int i = 0; i < 3; i++) { // check if they are numbers
                if (!tokens[i].empty() && tokens[i].find_first_not_of("0123456789-.Ee") != string::npos) {skip = true;}
            }
            if (skip) {continue;}

            // now we are most likely beyond any headers
            double _q, _I, _sigma;
            _q = std::stod(tokens[0]); // we know for sure that the strings are convertible to numbers (boost check)
            _I = std::stod(tokens[1]);
            _sigma = std::stod(tokens[2]);

            if (_q > 10) {continue;} // probably not a q-value if it's larger than 10

            // add the values to our vectors
            q.push_back(_q);
            I.push_back(_I);
            sigma.push_back(_sigma); 
        }
        return std::make_tuple(q, I, sigma);
    }
};