#include <utility/Dataset.h>
#include <math/SimpleLeastSquares.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

#include <vector>
#include <string>
#include <fstream>

using std::vector, std::string;

void IDataset::reduce(unsigned int target, bool log) {
    if (size() < target) {throw except::invalid_operation("Error in Dataset::reduce: Target cannot be larger than the size of the data set.");}
    Matrix<double> reduced(0, M);

    if (log) {
        double start = std::log10(x(0)); 
        double end = std::log10(x(size()-1));
        double width = (end - start)/target;

        unsigned int j = 0;
        for (unsigned int i = 0; i < size(); i++) {
            double val = std::log10(x(i));
            if (start + j*width < val) { // find the first x-value higher than our next sampling point
                reduced.push_back(row(i));
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
                reduced.push_back(row(i));
            }
        }
    }

    *this = std::move(reduced);
    options.draw_line = false;
    options.draw_markers = true;
}

void IDataset::limit_x(const Limit& limits) {
    if (size() == 0) {return;}
    if (limits.min < x(0) && x(size()-1) < limits.max) {return;}

    Matrix<double> limited(0, M); 
    for (unsigned int i = 0; i < size(); i++) {
        double val = x(i);
        if (val < limits.min || limits.max < val) {continue;}
        limited.push_back(row(i));
    }

    *this = std::move(limited);
}

void IDataset::limit_x(double min, double max) {limit_x({min, max});}

void IDataset::limit_y(const Limit& limits) {
    if (size() == 0) {return;}

    Matrix<double> limited(0, M);
    for (unsigned int i = 0; i < size(); i++) {
        double val = y(i);
        if (val < limits.min || limits.max < val) {continue;}
        limited.push_back(row(i));
    }

    *this = std::move(limited);
}

void IDataset::limit_y(double min, double max) {limit_y({min, max});}

void IDataset::operator=(const Matrix<double>&& other) {
    if (other.M != M) {throw except::invalid_operation("Error in Dataset::operator=: Matrix has wrong number of columns.");}
    this->data = std::move(other.data);
    this->N = other.N;
}

std::size_t IDataset::size() const noexcept {
    return N;
}

bool IDataset::is_logarithmic() const noexcept {
    // generate a new dataset containing exp(Deltax) and fit it with linear regression.
    // if the fit is decent, the data must have been logaritmic
    SimpleDataset exp_data;
    for (unsigned int i = 1; i < size(); i++) {
        exp_data.push_back(x(i), std::exp(x(i)-x(i-1)), 1);
    }

    SimpleLeastSquares fit(std::move(exp_data));
    auto res = fit.fit();

    std::cout << "DATASET IS_LOGARITHMIC FIT CHI: " << res->fval/res->dof << std::endl;
    std::cout << "chi: " << res->fval << std::endl;
    std::cout << "dof: " << res->dof << std::endl;
    return res->fval/res->dof < 10;
}

void IDataset::save(std::string path) const {
    utility::create_directories(path);

    // check if file was succesfully opened
    std::ofstream output(path);
    if (!output.is_open()) {throw std::ios_base::failure("Error in IntensityFitter::save: Could not open file \"" + path + "\"");}

    // prepare header & writer function
    std::function<string(unsigned int)> writer;
    string header = "broken header\n";
    if (M == 2) {
        header = "x y\n";
        writer = [this] (unsigned int i) {return std::to_string(index(i, 0)) + " " + std::to_string(index(i, 1)) + "\n";};
    } else if (M == 3) {
        header = "x y yerr\n";
        writer = [this] (unsigned int i) {return std::to_string(index(i, 0)) + " " + std::to_string(index(i, 1)) + " " + std::to_string(index(i, 2)) + "\n";};
    } else if (M == 4) {
        header = "x y yerr xerr\n";
        writer = [this] (unsigned int i) {return std::to_string(index(i, 0)) + " " + std::to_string(index(i, 1)) + " " + std::to_string(index(i, 2)) + " " + std::to_string(index(i, 3)) + "\n";};
    } else {
        throw except::invalid_operation("Error in IDataset::save: Dataset has wrong number of columns.");
    }

    // write to disk
    output << header;
    for (unsigned int i = 0; i < size(); i++) {
        output << writer(i);
    }
    output.close();
}

Limit IDataset::span_x() const noexcept {
    auto x = this->x();
    auto[min, max] = std::minmax_element(x.begin(), x.end());
    return Limit(*min, *max);
}

Limit IDataset::span_y() const noexcept {
    auto y = this->y();
    auto[min, max] = std::minmax_element(y.begin(), y.end());
    return Limit(*min, *max);
}

Limit IDataset::span_y_positive() const noexcept {
    auto y = this->y();
    if (size() == 0) {
        return Limit(0, 0);
    }

    Limit limits(y[0], y[0]);
    for (double val : y) {
        if (0 < val) {
            limits.min = std::min(val, limits.min);
        }
        limits.max = std::max(val, limits.max);
    }
    return limits;
}