#include <vector>

#include <Dataset.h>
#include <Exceptions.h>
#include <cmath>

using std::vector;

Dataset::Dataset(const std::vector<double>& x, const std::vector<double>& y) : x(x), y(y) {
    check_sizes();
}

Dataset::Dataset(const std::vector<double>& x, const std::vector<double>& y, const string xlabel, const string ylabel) 
    : xlabel(xlabel), ylabel(ylabel), x(x), y(y) {
        check_sizes();
}

void Dataset::reduce(unsigned int target) {
    if (size() < target) {throw except::invalid_operation("Error in Dataset::reduce: Target cannot be larger than the size of the data set.");}
    vector<double> new_x; new_x.reserve(target);
    vector<double> new_y; new_y.reserve(target);

    int ratio = std::floor(size()/target);
    for (unsigned int i = 0; i < size(); i++) {
        if (i % ratio == 0) {
            new_x.push_back(x[i]);
            new_y.push_back(y[i]);
        }
    }

    x = std::move(new_x);
    y = std::move(new_y);
}

std::size_t Dataset::size() const {return x.size();}

std::vector<double>& Dataset::get(const string label) {
    if (xlabel == label) {return x;}
    if (ylabel == label) {return y;}
    throw except::unknown_argument("Error in Dataset::get: No data is labelled as \"" + label + "\". Labels: " + xlabel + ", y: " + ylabel + ".");
}

void Dataset::check_sizes() const {
    if (x.size() != y.size()) {throw except::size_error("Error in Dataset::Dataset: x and y must have same size!");}
}