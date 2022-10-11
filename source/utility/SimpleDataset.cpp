#include <math/SimpleLeastSquares.h>
#include <math/Statistics.h>
#include <utility/SimpleDataset.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <utility/Settings.h>

#include <vector>
#include <string>
#include <fstream>
#include <random>

using std::vector, std::string;

SimpleDataset::SimpleDataset(unsigned int N, unsigned int M) : Dataset(N, M) {}

SimpleDataset::SimpleDataset(unsigned int rows) noexcept : Dataset(rows, 3) {}

SimpleDataset::SimpleDataset() noexcept : SimpleDataset(0) {}

SimpleDataset::SimpleDataset(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) : SimpleDataset(x.size()) {
    if (x.size() != y.size() || x.size() != yerr.size()) {
        throw except::size_error("SimpleDataset::SimpleDataset: x, y, and yerr must have the same size.");
    }
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], yerr[i]};
    }
}

SimpleDataset::SimpleDataset(std::vector<double> x, std::vector<double> y) : SimpleDataset(x, y, std::vector<double>(x.size())) {}

SimpleDataset::SimpleDataset(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel) : Dataset({x, y, std::vector<double>(x.size())}) {
    set_col_names({xlabel, ylabel, std::string(ylabel)+"err"});
    options.xlabel = xlabel;
    options.ylabel = ylabel;
}

// load is responsible for preparing the class
SimpleDataset::SimpleDataset(std::string path) {
    load(path);
}

void SimpleDataset::reduce(unsigned int target, bool log) {
    if (size() < target) {throw except::invalid_operation("Dataset::reduce: Target cannot be larger than the size of the data set. (" + std::to_string(target) + " > " + std::to_string(size()) + ")");}
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

void SimpleDataset::limit_x(const Limit& limits) {
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

void SimpleDataset::limit_x(double min, double max) {limit_x({min, max});}

void SimpleDataset::limit_y(const Limit& limits) {
    if (size() == 0) {return;}

    Matrix<double> limited(0, M);
    for (unsigned int i = 0; i < size(); i++) {
        double val = y(i);
        if (val < limits.min || limits.max < val) {continue;}
        limited.push_back(row(i));
    }

    *this = std::move(limited);
}

void SimpleDataset::limit_y(double min, double max) {limit_y({min, max});}

void SimpleDataset::operator=(const Matrix<double>&& other) {
    if (other.M != M) {throw except::invalid_operation("Dataset::operator=: Matrix has wrong number of columns.");}
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

    vector<double> x(size), y(size), yerr(size);
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
        utility::print_warning("Warning in SimpleDataset::simulate_errors: Dataset is empty.");
        return;
    }
    double y0 = y(0);
    auto x = this->x();
    auto y = this->y();
    auto yerr = this->yerr();
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return std::pow(y*x, 0.85);});
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return std::pow(y, 0.15)*std::pow(y0, 0.35)*std::pow(x, -0.85)/10000 + std::pow(x, 5)/100;});
    // std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return y/x*1e-4 + 1e-4;});
    std::transform(y.begin(), y.end(), x.begin(), yerr.begin(), [&y0] (double y, double x) {return y0/std::pow(x, 1.2)*1e-5 + 1e-4*y0;});    
}

Point2D SimpleDataset::get_point(unsigned int index) const {
    return Point2D(x(index), y(index), yerr(index));
}

Point2D SimpleDataset::find_minimum() const {
    if (size() == 0) {
        utility::print_warning("Warning in SimpleDataset::find_minimum: Dataset is empty.");
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
    SimpleDataset data; // rebinned dataset

    for (unsigned int i = 0; i < size(); i++) {
        // determine how many data points to fold into one
        unsigned int fold;
        if (0.1 < x(i)) {fold = 8;}
        else if (0.06 < x(i)) {fold = 4;}
        else if (0.03 < x(i)) {fold = 2;}
        else {fold = 1;}

        // loop over each data point to be folded
        double siginv = 0, sumw = 0, qsum = 0;
        unsigned int ss = 0;
        for (; (ss < fold) && (i < size()); ss++) {
            siginv += (std::pow(yerr(i), -2));
            sumw += y(i)/(std::pow(yerr(i), 2));
            qsum += x(i);
            i++;
        }

        // average their values into a single new one
        double q = qsum/ss;
        double I = sumw/siginv;
        double Ierr = std::pow(siginv, -0.5);
        data.push_back(q, I, Ierr);
    }
    *this = data;
}

void SimpleDataset::load(std::string path) {
    Dataset::load(path);
    if (M < 3) {
        throw except::io_error("SimpleDataset::load: Dataset has too few columns.");
    }
    else if (M > 3) {
        utility::print_warning("Warning in SimpleDataset::load: Dataset has " + std::to_string(M) + " columns, while this class only supports operations on 3. Ensure that the file is of the format [x | y | yerr | rest].");
    }
    names = {"q", "I", "Ierr", "qerr"}; // set column names
    unsigned int N = size();
    limit_x(setting::axes::qmin, setting::axes::qmax);
    if (N != size() && setting::general::verbose) {
        std::cout << "\tRemoved " << N - size() << " data points outside specified q-range [" << setting::axes::qmin << ", " << setting::axes::qmax << "]." << std::endl;
    }
}

void SimpleDataset::remove_consecutive_duplicates() {
    if (size() == 0) {
        utility::print_warning("Warning in SimpleDataset::remove_consecutive_duplicates: Dataset is empty.");
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

#include <math/MovingAverager.h>
void SimpleDataset::moving_average(unsigned int window_size) {
    auto newy = MovingAverage::average_half(y(), window_size);

    unsigned int offset = (window_size - 1)/2;
    Matrix<double> newdata(N-window_size+1, M);
    for (unsigned int i = 0; i < newdata.N; i++) {
        auto row = this->row(i+offset);
        row[1] = newy[i];
        newdata.row(i) = this->row(i + offset);
    }

    this->data = newdata.data;
    this->N = newdata.N;
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