#include <vector>
#include <cmath>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include <TGraph.h>
#include <TGraphErrors.h>

#include <Dataset.h>
#include <Exceptions.h>
#include <settings.h>

using std::vector, std::string;

Dataset::Dataset(const string file) {
    read(file);
}

Dataset::Dataset(const std::vector<double>& x, const std::vector<double>& y) : x(x), y(y) {
    validate_sizes();
}

Dataset::Dataset(const std::vector<double>& x, const std::vector<double>& y, const string xlabel, const string ylabel) 
    : xlabel(xlabel), ylabel(ylabel), xerrlabel(xlabel+"err"), yerrlabel(ylabel+"err"), x(x), y(y) {
        validate_sizes();
}

Dataset::Dataset(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& yerr) : x(x), y(y), yerr(yerr) {
    xerr.resize(yerr.size());
    validate_sizes();
}

Dataset::Dataset(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xerr, const std::vector<double>& yerr)
    : x(x), y(y), xerr(xerr), yerr(yerr) {
        validate_sizes();
    }

void Dataset::reduce(unsigned int target, bool log) {
    if (size() < target) {throw except::invalid_operation("Error in Dataset::reduce: Target cannot be larger than the size of the data set.");}
    vector<double> new_x; new_x.reserve(target);
    vector<double> new_y; new_y.reserve(target);

    if (log) {
        double start = std::log10(x[0]); 
        double end = std::log10(x[x.size()-1]);
        double width = (end - start)/target;

        unsigned int j = 0;
        for (unsigned int i = 0; i < size(); i++) {
            double val = std::log10(x[i]);
            if (start + j*width < val) { // find the first x-value higher than our next sampling point
                new_x.push_back(x[i]);
                new_y.push_back(y[i]);
                j++;
            }
            while (start + j*width < val) { // it may be necessary to skip a few sampled points, especially at the beginning
                j++;
            }
        }
    } else {
        int ratio = std::floor(size()/target);
        for (unsigned int i = 0; i < size(); i++) {
            if (i % ratio == 0) {
                new_x.push_back(x[i]);
                new_y.push_back(y[i]);
            }
        }
    }

    x = std::move(new_x);
    y = std::move(new_y);

    draw_as_line = false;
}

std::size_t Dataset::size() const {return x.size();}

void Dataset::limit(const Limit& limits) {
    if (limits.min < x[0] && x[size()-1] < limits.max) {return;}

    vector<double> new_x; new_x.reserve(size());
    vector<double> new_y; new_y.reserve(size());

    for (unsigned int i = 0; i < size(); i++) {
        double val = x[i];
        if (val < limits.min || limits.max < val) {continue;}
        new_x.push_back(x[i]);
        new_y.push_back(y[i]);
    }

    x = std::move(new_x);
    y = std::move(new_y);
}

std::vector<double>& Dataset::get(const string label) {
    if (xlabel == label) {return x;}
    if (ylabel == label) {return y;}
    if (xerrlabel == label) {return xerr;}
    if (yerrlabel == label) {return yerr;}
    throw except::unknown_argument("Error in Dataset::get: No data is labelled as \"" + label + "\". Labels: " + xlabel + ", y: " + ylabel + ".");
}

void Dataset::validate_sizes() const {
    if (x.size() != y.size()
        || (!yerr.empty() && yerr.size() != x.size())
        || (!xerr.empty() && xerr.size() != x.size())) 
        {
            throw except::size_error("Error in Dataset::Dataset: x and y must have same size!");
    }
}

std::unique_ptr<TGraph> Dataset::plot() const {
    if (!yerr.empty()) {return std::make_unique<TGraphErrors>(x.size(), x.data(), y.data(), xerr.data(), yerr.data());}
    else {return std::make_unique<TGraph>(x.size(), x.data(), y.data());}
}

void Dataset::scale_errors(double factor) {
    if (!xerr.empty()) {std::transform(xerr.begin(), xerr.end(), xerr.begin(), [&factor] (double val) {return factor*val;});}
    if (!yerr.empty()) {std::transform(yerr.begin(), yerr.end(), yerr.begin(), [&factor] (double val) {return factor*val;});}
}

void Dataset::scale_y(double factor) {
    std::transform(y.begin(), y.end(), y.begin(), [&factor] (double val) {return val*factor;});
    if (!yerr.empty()) {std::transform(yerr.begin(), yerr.end(), yerr.begin(), [&factor] (double val) {return factor*val;});}
}

void SAXSDataset::simulate_errors() {
    if (yerr.empty()) {yerr.resize(size());}
    if (xerr.empty()) {xerr.resize(size());}

    double base = y[0]+1;
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&base] (double y, double x) {return std::pow(y, 0.75)*x;});
    std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&base] (double y, double x) {return std::pow(y*x, 0.85);});
}

void SAXSDataset::set_resolution(unsigned int resolution) {
    this->resolution = resolution;
    limit(Limit(0, 2*M_PI/resolution));
}

void Dataset::read(const string file) {
    // check if file was succesfully opened
    std::ifstream input(file);
    if (!input.is_open()) {throw std::ios_base::failure("Error in IntensityFitter::read: Could not open file \"" + file + "\"");}

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

        // check user-defined limits
        if (_q < setting::fit::q_low) {continue;}
        if (_q > setting::fit::q_high) {continue;}

        // add the values to our vectors
        x.push_back(_q);
        y.push_back(_I);
        yerr.push_back(_sigma); 
    }
}