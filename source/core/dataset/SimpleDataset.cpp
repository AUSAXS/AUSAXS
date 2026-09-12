// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <dataset/SimpleDataset.h>

#include <dataset/DatasetFactory.h>
#include <hist/Histogram.h>
#include <math/Statistics.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>
#include <utility/Random.h>

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

#ifndef NDEBUG
    #include <iostream>  // only the asserts below print
#endif

using namespace ausaxs;

SimpleDataset::SimpleDataset() : SimpleDataset(0) {}
SimpleDataset::SimpleDataset(const SimpleDataset& d) = default;
SimpleDataset::SimpleDataset(SimpleDataset&& d) noexcept = default;
SimpleDataset& SimpleDataset::operator=(const SimpleDataset& other) = default;
SimpleDataset& SimpleDataset::operator=(SimpleDataset&& other) noexcept = default;
SimpleDataset::~SimpleDataset() = default;

SimpleDataset::SimpleDataset(const Dataset& d) : SimpleDataset(d.size()) {
    if (d.data.M <= 1) {throw except::invalid_argument("SimpleDataset::SimpleDataset: Dataset must have at least two columns.");}

    if (d.data.M == 3) {
        data = d.data;
    } else {
        for (int i = 0; i < data.N; i++) {
            row(i) = {d.x(i), d.y(i), 0};
        }
    }
}

SimpleDataset::SimpleDataset(const hist::Histogram& h) : SimpleDataset(h.as_dataset()) {}

SimpleDataset::SimpleDataset(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& yerr) : SimpleDataset(static_cast<int>(x.size())) {initialize(x, y, yerr);}

SimpleDataset::SimpleDataset(const std::vector<double>& x, const std::vector<double>& y) : SimpleDataset(static_cast<int>(x.size())) {initialize(x, y);}

SimpleDataset::SimpleDataset(int N, int M) : Dataset(N, M) {}

SimpleDataset::SimpleDataset(int rows) noexcept : Dataset(rows, 3) {}

SimpleDataset::SimpleDataset(const io::ExistingFile& path) : SimpleDataset() {
    auto data = factory::DatasetFactory::construct(path, 3);
    *this = std::move(*data);
}

void SimpleDataset::initialize(const std::vector<double>& x, const std::vector<double>& y) {
    assert([&]() -> bool {
        if (x.size() == y.size()) {return true;}
        std::cout << "SimpleDataset::initialize: x and y must have the same size (" << x.size() << ", " << y.size() << ")." << std::endl;
        return false;
    }() && "SimpleDataset::initialize: x and y must have the same size.");
    for (int i = 0; i < static_cast<int>(x.size()); i++) {
        row(i) = {x[i], y[i], 0};
    }
}

void SimpleDataset::initialize(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& yerr) {
    assert([&]() -> bool {
        if (x.size() == y.size() && x.size() == yerr.size()) {return true;}
        std::cout << "SimpleDataset::initialize: x, y, and yerr must have the same size (" << x.size() << ", " << y.size() << ", " << yerr.size() << ")." << std::endl;
        return false;
    }() && "SimpleDataset::initialize: x, y, and yerr must have the same size.");
    for (int i = 0; i < static_cast<int>(x.size()); i++) {
        row(i) = {x[i], y[i], yerr[i]};
    }
}

void SimpleDataset::reduce(int target, bool log) {
    if (size() < target) {
        if (settings::general::verbose) {
            console::print_warning("Warning in SimpleDataset::reduce: Dataset is already smaller than target size.");
        } 
        return;
    }

    Matrix<double> reduced(0, data.M);

    if (log) {
        double start = std::log10(x().front()); 
        double end = std::log10(x().back());
        double width = (end - start)/(target-1);

        reduced.push_back(row(0));
        int j = 1;
        for (int i = 1; i < size(); i++) {
            double val = std::log10(x(i));

            // find the first x-value higher than our next sampling point
            if (start + j*width <= val) { 
                reduced.push_back(row(i));
                j++;
            }
            
            // it may be necessary to skip a few sampled points, especially at the beginning            
            while (start + j*width < val) { 
                j++;
            }
        }
    } else {
        int j = 0;
        double ratio = double(size())/target;
        for (int i = 0; i < size(); i++) {
            if (i >= j*ratio) {
                reduced.push_back(row(i));
                j++;
            }
        }
    }
    *this = std::move(reduced);
}

SimpleDataset& SimpleDataset::operator=(Matrix<double>&& other) { // NOLINT - only the matrix contents are moved
    if (other.M != data.M) {throw except::invalid_operation("SimpleDataset::operator=: Matrix has wrong number of columns.");}
    this->data.data = std::move(other.data);
    this->data.N = other.N;
    return *this;
}

Limit SimpleDataset::span_x() const noexcept {
    if (empty()) {
        return {0, 0};
    }
    auto x = this->x();
    auto[min, max] = std::ranges::minmax_element(x);
    return {*min, *max};
}

Limit SimpleDataset::span_y() const noexcept {
    if (empty()) {
        return {0, 0};
    }
    auto y = this->y();
    auto[min, max] = std::ranges::minmax_element(y);
    return {*min, *max};
}

Limit SimpleDataset::get_xlimits() const noexcept {return span_x();}

Limit SimpleDataset::get_ylimits() const noexcept {return span_y();}

Limit SimpleDataset::span_y_positive() const noexcept {
    auto y = this->y();
    if (empty()) {
        return {0, 0};
    }

    Limit limits;
    // find first non-zero y value
    int i = 0;
    for (; i < size(); i++) {
        if (0 < y[i]) {
            limits.min = y[i];
            limits.max = y[i];
            break;
        }
    }

    // continue search for lower mins and higher max
    for (; i < size(); i++) {
        double val = y[i];
        if (0 < val) {
            limits.min = std::min(val, limits.min);
        }
        limits.max = std::max(val, limits.max);
    }
    return limits;
}

SimpleDataset SimpleDataset::generate_random_data(int size, double val) {
    return generate_random_data(size, -val, val);
}

SimpleDataset SimpleDataset::generate_random_data(int size, double min, double max) {
    auto uniform = std::uniform_real_distribution<double>(min, max);

    std::vector<double> x(size), y(size), yerr(size);
    for (int i = 0; i < size; i++) {
        x[i] = i;
        y[i] = uniform(random::generator());
        yerr[i] = y[i]*0.1;
    }
    return {x, y, yerr};
}

void SimpleDataset::push_back(double x, double y, double yerr) {
    data.extend(1);
    row(data.N-1) = {x, y, yerr};
}

double SimpleDataset::normalize(double y0) {
    double scale = y0/y(0);
    scale_y(scale);
    return scale;
}

void SimpleDataset::scale_errors(double factor) {
    auto yerr = this->yerr();
    std::ranges::transform(yerr, yerr.begin(), [&factor] (double val) {return factor*val;});
}

void SimpleDataset::scale_y(double factor) {
    auto y = this->y();
    auto yerr = this->yerr();
    std::ranges::transform(y, y.begin(), [&factor] (double val) {return val*factor;});
    std::ranges::transform(yerr, yerr.begin(), [&factor] (double val) {return factor*val;});
}

void SimpleDataset::simulate_noise() {
    auto fun = [] (double y, double yerr) {
        return random::gaussian(y, yerr);
    };

    auto y = this->y();
    auto yerr = this->yerr();
    std::ranges::transform(y, yerr, y.begin(), fun);
}

void SimpleDataset::simulate_errors() {
    if (empty()) {
        if (settings::general::verbose) {
            console::print_warning("Warning in SimpleDataset::simulate_errors: Dataset is empty.");
        }
        return;
    }
    double y0 = y(0);
    auto x = this->x();
    auto y = this->y();
    auto yerr = this->yerr();
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return std::pow(y*x, 0.85);});
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return std::pow(y, 0.15)*std::pow(y0, 0.35)*std::pow(x, -0.85)/10000 + std::pow(x, 5)/100;});
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return y/x*1e-4 + 1e-4;});
    std::ranges::transform(y, x, yerr.begin(), [&y0] (double, double x) {return y0*(1 + 0.1/std::pow(x, 1.2))*1e-4;});
}

Point2D SimpleDataset::get_point(int index) const {
    if (data.M < 3) {return {x(index), y(index)};}
    return {x(index), y(index), yerr(index)};
}

Point2D SimpleDataset::find_minimum() const {
    auto res = Dataset::find_minimum(1);
    return {res[0], res[1], res[2]};
}

void SimpleDataset::push_back(const Point2D& point) noexcept {
    push_back(point.x, point.y, point.yerr);
}

void SimpleDataset::rebin() noexcept {
    SimpleDataset newdata; // rebinned dataset

    std::function<void(int, int&)> func;
    if (std::accumulate(yerr().begin(), yerr().end(), 0.0) == 0) {
        func = [&newdata, this] (int nfold, int& index) {
            double wsum = 0, qsum = 0, folds = 0;
            for (; (folds < nfold) && (index < size()); folds++) {
                wsum += y(index);
                qsum += x(index++);
            }
            assert(folds != 0 && "SimpleDataset::rebin: Division by zero.");
            newdata.push_back(qsum/folds, wsum/folds, 0);
        };
    } else {
        func = [&newdata, this] (int nfold, int& index) {
            double siginv = 0, wsum = 0, qsum = 0, folds = 0;
            for (; (folds < nfold) && (index < size()); folds++) {
                siginv += (std::pow(yerr(index), -2));
                wsum += y(index)/(std::pow(yerr(index), 2));
                qsum += x(index++);
            }
            assert(siginv != 0 && "SimpleDataset::rebin: Division by zero.");
            newdata.push_back(qsum/folds, wsum/siginv, std::pow(siginv, -0.5));
        };
    }

    // note: func() advances 'i' past every point it consumed, so the loop must not increment it again
    for (int i = 0; i < size(); ) {
        // determine how many data points to fold into one
        int fold;
        if (0.1 < x(i)) {fold = 8;}
        else if (0.06 < x(i)) {fold = 4;}
        else if (0.03 < x(i)) {fold = 2;}
        else {fold = 1;}

        // fold data points
        func(fold, i);
    }
    console::print_text("Rebinned dataset from " + std::to_string(size()) + " to " + std::to_string(newdata.size()) + " data points.");
    *this = std::move(newdata);
}

void SimpleDataset::load(const io::ExistingFile& path) {
    Dataset::load(path);
}

void SimpleDataset::remove_consecutive_duplicates() {
    if (empty()) {
        console::print_warning("Warning in SimpleDataset::remove_consecutive_duplicates: Dataset is empty.");
        return;
    }

    Matrix<double> new_data(data.N, data.M);
    new_data.row(0) = this->row(0);

    int index = 1;
    double v = y(0);
    for (int i = 1; i < size(); i++) {
        if (y(i) != v) {
            new_data.row(index++) = this->row(i);
            v = y(i);
        }
    }
    new_data.resize(index, data.M);
    this->assign_matrix(std::move(new_data));
}

double SimpleDataset::mean() const {
    return stats::mean(y());
}

double SimpleDataset::weighted_mean() const {
    return stats::weighted_mean(y(), yerr());
}

double SimpleDataset::std() const {
    return stats::std(y());
}

double SimpleDataset::weighted_mean_error() const {
    return stats::weighted_mean_error(yerr());
}

// bool SimpleDataset::operator==(const SimpleDataset& other) const = default;

bool SimpleDataset::operator==(const SimpleDataset& other) const {
    if (size() != other.size()) {return false;}
    if (this->data != other.data) {return false;}
    return true;
}