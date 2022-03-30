#include <vector>

#include <Dataset.h>
#include <Exceptions.h>
#include <cmath>

#include <TGraph.h>
#include <TGraphErrors.h>

using std::vector;

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
}

std::size_t Dataset::size() const {return x.size();}

void Dataset::limit(const Limit& limits) {
    if (limits.min < x[0] && x[size()-1] < limits.max) {return;}

    vector<double> new_x; new_x.reserve(size());
    vector<double> new_y; new_y.reserve(size());

    unsigned int i = 0;
    while (x[i] <= limits.max) {
        new_x.push_back(x[i]);
        new_y.push_back(y[i++]);
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

void SAXSDataset::simulate_errors() {
    if (yerr.empty()) {yerr.resize(size());}
    if (xerr.empty()) {xerr.resize(size());}

    double base = y[0]+1;
    // std::transform(y.begin(), y.end(), yerr.begin(), [&base] (double val) {return log10(base-val) + 1;});
    // std::transform(y.begin(), y.end(), yerr.begin(), [&base] (double val) {return 20*log10(base-val)*std::sqrt(val) + 1;});
    std::transform(y.begin(), y.end(), yerr.begin(), [&base] (double val) {return std::sqrt(val);});
}

void SAXSDataset::set_resolution(unsigned int resolution) {
    this->resolution = resolution;
    limit(Limit(0, 2*M_PI/resolution));
}