#include <math/SimpleLeastSquares.h>
#include <math/Statistics.h>
#include <dataset/SimpleDataset.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <dataset/DatasetFactory.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>

#include <vector>
#include <string>
#include <fstream>
#include <random>

SimpleDataset::SimpleDataset(const Dataset& d) : SimpleDataset(d.size()) {
    if (d.M <= 1) {
        throw except::invalid_argument("SimpleDataset::SimpleDataset: Dataset must have at least two columns.");
    } else if (d.M == 3) {
        data = d.data;
    } else {
        for (unsigned int i = 0; i < N; i++) {
            row(i) = {d.x(i), d.y(i), 0};
        }
    }
}

SimpleDataset::SimpleDataset(unsigned int N, unsigned int M) : Dataset(N, M) {}

SimpleDataset::SimpleDataset(unsigned int rows) noexcept : Dataset(rows, 3) {}

SimpleDataset::SimpleDataset() noexcept : SimpleDataset(0) {}

SimpleDataset::SimpleDataset(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) : SimpleDataset(x.size()) {
    if (x.size() != y.size() || x.size() != yerr.size()) {
        throw except::size_error("SimpleDataset::SimpleDataset: x, y, and yerr must have the same size (" + std::to_string(x.size()) + ", " + std::to_string(y.size()) + ", " + std::to_string(yerr.size()) + ".");
    }
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], yerr[i]};
    }
}

SimpleDataset::SimpleDataset(std::vector<double> x, std::vector<double> y) : SimpleDataset(x, y, std::vector<double>(y.size())) {}

SimpleDataset::SimpleDataset(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel) : Dataset({x, y, std::vector<double>(y.size())}) {
    set_col_names({xlabel, ylabel, std::string(ylabel)+"err"});
    options.xlabel = xlabel;
    options.ylabel = ylabel;
}

SimpleDataset::SimpleDataset(std::string path) : SimpleDataset() {
    auto data = factory::DatasetFactory::construct(path, 3);
    *this = std::move(*data);
}

void SimpleDataset::reduce(unsigned int target, bool log) {
    if (size() < target) {throw except::invalid_operation("SimpleDataset::reduce: Target cannot be larger than the size of the data set. (" + std::to_string(target) + " > " + std::to_string(size()) + ")");}
    Matrix<double> reduced(0, M);

    if (log) {
        double start = std::log10(x(0)); 
        double end = std::log10(x().back());
        double width = (end - start)/target;

        reduced.push_back(row(0));
        unsigned int j = 1;
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

void SimpleDataset::operator=(const Matrix<double>&& other) {
    if (other.M != M) {throw except::invalid_operation("SimpleDataset::operator=: Matrix has wrong number of columns.");}
    this->data = std::move(other.data);
    this->N = other.N;
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
    extend(1);
    row(N-1) = {x, y, yerr};
}

void SimpleDataset::normalize(double y0) {
    scale_y(y0/y(0));
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
        console::print_warning("Warning in SimpleDataset::simulate_errors: Dataset is empty.");
        return;
    }
    double y0 = y(0);
    auto x = this->x();
    auto y = this->y();
    auto yerr = this->yerr();
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return std::pow(y*x, 0.85);});
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return std::pow(y, 0.15)*std::pow(y0, 0.35)*std::pow(x, -0.85)/10000 + std::pow(x, 5)/100;});
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return y/x*1e-4 + 1e-4;});
    std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double, double x) {return y0/std::pow(x, 1.2)*1e-5 + 1e-4*y0;});    
}

Point2D SimpleDataset::get_point(unsigned int index) const {
    if (M < 3) {return Point2D(x(index), y(index));}
    else {      return Point2D(x(index), y(index), yerr(index));}
}

Point2D SimpleDataset::find_minimum() const {
    if (size() == 0) {
        console::print_warning("Warning in SimpleDataset::find_minimum: Dataset is empty.");
        return Point2D(0, 0, 0);
    }
    
    unsigned int min_index = 0;
    double min_value = y(0);
    for (unsigned int i = 1; i < size(); i++) {
        if (y(i) < min_value) {
            min_index = i;
            min_value = y(i);
        }
    }
    return get_point(min_index);
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

    Matrix new_data(N, M);
    new_data.row(0) = this->row(0);

    unsigned int index = 1;
    double v = y(0);
    for (unsigned int i = 1; i < size(); i++) {
        if (y(i) != v) {
            new_data.row(index++) = this->row(i);
            v = y(i);
        }
    }
    new_data.resize(index, M);
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

void SimpleDataset::sort_x() {
    Matrix<double> newdata(N, M);
    std::vector<unsigned int> indices(N);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [this] (unsigned int i, unsigned int j) {return x(i) < x(j);});
    for (unsigned int i = 0; i < N; i++) {
        newdata.row(i) = this->row(indices[i]);
    }
    this->assign_matrix(std::move(newdata));
}
