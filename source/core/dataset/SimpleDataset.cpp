// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <dataset/SimpleDataset.h>
#include <math/Statistics.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <dataset/DatasetFactory.h>
#include <hist/Histogram.h>
#include <settings/GeneralSettings.h>

#include <vector>
#include <string>
#include <random>

using namespace ausaxs;

SimpleDataset::SimpleDataset() : SimpleDataset(0) {}
SimpleDataset::SimpleDataset(const SimpleDataset& d) = default;
SimpleDataset::SimpleDataset(SimpleDataset&& d) = default;
SimpleDataset& SimpleDataset::operator=(const SimpleDataset& other) = default;
SimpleDataset& SimpleDataset::operator=(SimpleDataset&& other) noexcept = default;
SimpleDataset::~SimpleDataset() = default;

SimpleDataset::SimpleDataset(const Dataset& d) : SimpleDataset(d.size()) {
    if (d.data.M <= 1) {
        throw except::invalid_argument("SimpleDataset::SimpleDataset: Dataset must have at least two columns.");
    } else if (d.data.M == 3) {
        data = d.data;
    } else {
        for (unsigned int i = 0; i < data.N; i++) {
            row(i) = {d.x(i), d.y(i), 0};
        }
    }
}

SimpleDataset::SimpleDataset(const hist::Histogram& h) : SimpleDataset(h.as_dataset()) {}

SimpleDataset::SimpleDataset(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& yerr) : SimpleDataset(x.size()) {initialize(x, y, yerr);}

SimpleDataset::SimpleDataset(const std::vector<double>& x, const std::vector<double>& y) : SimpleDataset(x.size()) {initialize(x, y);}

SimpleDataset::SimpleDataset(const std::vector<double>& x, const std::vector<double>& y, std::string xlabel, std::string ylabel) : SimpleDataset(x.size()) {
    initialize(x, y);
    set_col_names({std::move(xlabel), ylabel, std::move(ylabel)+"err"});
}

SimpleDataset::SimpleDataset(unsigned int N, unsigned int M) : Dataset(N, M) {}

SimpleDataset::SimpleDataset(unsigned int rows) noexcept : Dataset(rows, 3) {}

SimpleDataset::SimpleDataset(const io::ExistingFile& path) : SimpleDataset() {
    auto data = factory::DatasetFactory::construct(path, 3);
    *this = std::move(*data);
}

void SimpleDataset::initialize(const std::vector<double>& x, const std::vector<double>& y) {
    #if DEBUG
        if (x.size() != y.size()) {
            throw except::size_error("SimpleDataset::SimpleDataset: x and y must have the same size (" + std::to_string(x.size()) + ", " + std::to_string(y.size()) + ").");
        }
    #endif
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], 0};
    }
}

void SimpleDataset::initialize(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& yerr) {
    #if DEBUG
        if (x.size() != y.size() || x.size() != yerr.size()) {
            throw except::size_error("SimpleDataset::SimpleDataset: x, y, and yerr must have the same size (" + std::to_string(x.size()) + ", " + std::to_string(y.size()) + ", " + std::to_string(yerr.size()) + ".");
        }
    #endif
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], yerr[i]};
    }
}

void SimpleDataset::reduce(unsigned int target, bool log) {
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
        unsigned int j = 1;
        for (unsigned int i = 1; i < size(); i++) {
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
        unsigned int j = 0;
        double ratio = double(size())/target;
        for (unsigned int i = 0; i < size(); i++) {
            if (i >= j*ratio) {
                reduced.push_back(row(i));
                j++;
            }
        }
    }
    *this = std::move(reduced);
}

void SimpleDataset::operator=(Matrix<double>&& other) {
    if (other.M != data.M) {throw except::invalid_operation("SimpleDataset::operator=: Matrix has wrong number of columns.");}
    this->data.data = std::move(other.data);
    this->data.N = other.N;
}

Limit SimpleDataset::span_x() const noexcept {
    if (size() == 0) {
        return Limit(0, 0);
    }
    auto x = this->x();
    auto[min, max] = std::minmax_element(x.begin(), x.end());
    return Limit(*min, *max);
}

Limit SimpleDataset::span_y() const noexcept {
    if (size() == 0) {
        return Limit(0, 0);
    }
    auto y = this->y();
    auto[min, max] = std::minmax_element(y.begin(), y.end());
    return Limit(*min, *max);
}

Limit SimpleDataset::get_xlimits() const noexcept {return span_x();}

Limit SimpleDataset::get_ylimits() const noexcept {return span_y();}

Limit SimpleDataset::span_y_positive() const noexcept {
    auto y = this->y();
    if (size() == 0) {
        return Limit(0, 0);
    }

    Limit limits;
    // find first non-zero y value
    unsigned int i = 0;
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

SimpleDataset SimpleDataset::generate_random_data(unsigned int size, double val) {
    return generate_random_data(size, -val, val);
}

SimpleDataset SimpleDataset::generate_random_data(unsigned int size, double min, double max) {
    std::random_device dev;
    std::mt19937 gen(dev());
    auto uniform = std::uniform_real_distribution<double>(min, max);

    std::vector<double> x(size), y(size), yerr(size);
    for (unsigned int i = 0; i < size; i++) {
        x[i] = i;
        y[i] = uniform(gen);
        yerr[i] = y[i]*0.1;
    }
    return SimpleDataset(x, y, yerr);
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
    std::transform(yerr.begin(), yerr.end(), yerr.begin(), [&factor] (double val) {return factor*val;});
}

void SimpleDataset::scale_y(double factor) {
    auto y = this->y();
    auto yerr = this->yerr();
    std::transform(y.begin(), y.end(), y.begin(), [&factor] (double val) {return val*factor;});
    std::transform(yerr.begin(), yerr.end(), yerr.begin(), [&factor] (double val) {return factor*val;});
}

void SimpleDataset::simulate_noise() {
    std::random_device dev;
    std::mt19937 gen(dev());
    auto fun = [&] (double y, double yerr) {
        auto gauss = std::normal_distribution<double>(y, yerr);
        return gauss(gen);
    };

    auto y = this->y();
    auto yerr = this->yerr();
    std::transform(y.begin(), y.end(), yerr.begin(), y.begin(), fun);
}

void SimpleDataset::simulate_errors() {
    if (size() == 0) {
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
    std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double, double x) {return y0*(1 + 0.1/std::pow(x, 1.2))*1e-4;});
}

Point2D SimpleDataset::get_point(unsigned int index) const {
    if (data.M < 3) {return Point2D(x(index), y(index));}
    else {           return Point2D(x(index), y(index), yerr(index));}
}

Point2D SimpleDataset::find_minimum() const {
    auto res = Dataset::find_minimum(1);
    return Point2D(res[0], res[1], res[2]);
}

void SimpleDataset::push_back(const Point2D& point) noexcept {
    push_back(point.x, point.y, point.yerr);
}

void SimpleDataset::rebin() noexcept {
    SimpleDataset newdata; // rebinned dataset

    std::function<void(unsigned int, unsigned int&)> func;
    if (std::accumulate(yerr().begin(), yerr().end(), 0.0) == 0) {
        func = [&newdata, this] (unsigned int nfold, unsigned int& index) {
            double wsum = 0, qsum = 0, folds = 0;
            for (; (folds < nfold) && (index < size()); folds++) {
                wsum += y(index);
                qsum += x(index++);
            }
            assert(folds != 0 && "SimpleDataset::rebin: Division by zero.");
            newdata.push_back(qsum/folds, wsum/folds, 0);
        };
    } else {
        func = [&newdata, this] (unsigned int nfold, unsigned int& index) {
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

    for (unsigned int i = 0; i < size(); i++) {
        // determine how many data points to fold into one
        unsigned int fold;
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
    names = {"q", "I", "Ierr", "qerr"}; // set column names
}

void SimpleDataset::remove_consecutive_duplicates() {
    if (size() == 0) {
        console::print_warning("Warning in SimpleDataset::remove_consecutive_duplicates: Dataset is empty.");
        return;
    }

    Matrix<double> new_data(data.N, data.M);
    new_data.row(0) = this->row(0);

    unsigned int index = 1;
    double v = y(0);
    for (unsigned int i = 1; i < size(); i++) {
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