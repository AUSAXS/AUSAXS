#include <math/SimpleLeastSquares.h>
#include <utility/Dataset.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <utility/Settings.h>

#include <vector>
#include <string>
#include <fstream>
#include <random>

using std::vector, std::string;

Column<double> IDataset::get(std::string column) {
    if (column == options.xlabel) {return x();}
    else if (column == options.ylabel) {return y();}
    else {throw except::invalid_argument("Error in IDataset::get: Column name \"" + column + "\" not recognized.");}
}

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
        exp_data.push_back(x(i), std::exp(x(i)-x(i-1)));
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

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
void IDataset::load(std::string path) {
    // check if file was succesfully opened
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("Error in IDataset::load: Could not open file \"" + path + "\"");}

    string line; // placeholder for the current line
    while(getline(input, line)) {
        if (line[0] == ' ') {line = line.substr(1);} // fix leading space
        vector<string> tokens;
        boost::split(tokens, line, boost::is_any_of(" ,\t")); // spaces, commas, and tabs can all be used as separators (but not a mix of them)

        // remove empty tokens
        for (unsigned int i = 0; i < tokens.size(); i++) {
            if (tokens[i].empty()) {tokens.erase(tokens.begin()+i);}
        }

        // determine if we are in some sort of header
        if (tokens.size() < 2 || tokens.size() > 4) {continue;} // too many separators
        bool skip = false;
        for (unsigned int i = 0; i < tokens.size(); i++) { // check if all tokens are numbers
            if (tokens[i].find_first_not_of("0123456789-.Ee\n\r") != string::npos) {
                skip = true;
            }
        }
        if (skip) {continue;}

        // now we are most likely beyond any headers
        double _q, _I, _sigma;
        _q = std::stod(tokens[0]); // we know for sure that the strings are convertible to numbers (boost check)
        _I = std::stod(tokens[1]);
        if (_q > 10) {continue;} // probably not a q-value if it's larger than 10

        // check user-defined limits
        if (_q < setting::fit::q_low) {continue;}
        if (_q > setting::fit::q_high) {continue;}

        // add the values to our vectors
        // this is a fair bit more complicated than strictly necessary
        if (tokens.size() == 4) {
            if (M == 4) {
                push_back({_q, _I, std::stod(tokens[2]), std::stod(tokens[3])});
            } else if (M == 3) {
                throw except::invalid_operation("Error in IDataset::load: File has four columns, but a SimpleDataset only supports three. Use a Dataset instance instead.");
            } else {
                throw except::unexpected("Error in IDataset::load: Unknown data layout.");
            }
        }
        else if (tokens.size() == 3) {
            if (M == 4) {
                push_back({_q, _I, std::stod(tokens[2]), 0});
            } else if (M == 3) {
                push_back({_q, _I, std::stod(tokens[2])});
            } else {
                throw except::unexpected("Error in IDataset::load: Unknown data layout.");
            }
        }
        else if (tokens.size() == 2) {
            if (M == 4) {
                push_back({_q, _I, 0, 0});
            } else if (M == 3) {
                push_back({_q, _I, 0});
            } else {
                throw except::unexpected("Error in IDataset::load: Unknown data layout.");
            }
        } else {
            throw except::unexpected("Error in IDataset::load: Unknown data layout.");
        }
    }
    input.close();
    if (size() == 0) {throw except::unexpected("Error in IDataset::load: No data was read from the file.");}
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

void SimpleDataset::name_columns(std::string xlabel, std::string ylabel) {
    options.xlabel = xlabel;
    options.ylabel = ylabel;
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

Point2D SimpleDataset::get_point(unsigned int index) const noexcept {
    return Point2D(x(index), y(index), yerr(index));
}

Point2D SimpleDataset::find_minimum() const noexcept {
    if (size() == 0) {
        utility::print_warning("Warning in SimpleDataset::find_minimum: Dataset is empty.");
        return Point2D(0, 0, 0);
    }
    auto y = this->y();
    auto it = std::min_element(y.begin(), y.end());
    unsigned int index = it - y.begin();
    return get_point(index);
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

        std::cout << "now folding " << i << " to " << i + fold << std::endl;

        // loop over each data point to be folded
        double siginv = 0, sumw = 0, qsum = 0;
        unsigned int ss = 0;
        for (; ss < fold; ss++) {
            std::cout << "checkpoint1" << std::endl;
            if (i == size()) {break;}
            std::cout << "checkpoint1" << std::endl;
            siginv += (std::pow(yerr(i), -2));
            std::cout << "checkpoint1" << std::endl;
            sumw += y(i)/(std::pow(yerr(i), 2));
            std::cout << "checkpoint1" << std::endl;
            qsum += x(i);
            std::cout << "checkpoint1" << std::endl;
            i++;
        }

        // average their values into a single new one
        double q = qsum/ss;
        double I = sumw/siginv;
        double Ierr = std::pow(siginv, -0.5);
        data.push_back(Point2D(q, I, Ierr));
    }
    data.save("temp/dataset/test.dat");
    *this = data;
}

void Dataset::scale_errors(double factor) {
    auto xerr = this->xerr();
    auto yerr = this->yerr();
    std::transform(xerr.begin(), xerr.end(), xerr.begin(), [&factor] (double val) {return factor*val;});
    std::transform(yerr.begin(), yerr.end(), yerr.begin(), [&factor] (double val) {return factor*val;});
}

void Dataset::push_back(double x, double y, double xerr, double yerr) {
    extend(1);
    row(N-1) = {x, y, yerr, xerr};
}

void Dataset::push_back(double x, double y) {
    push_back(x, y, 0, 0);
}

void Dataset::push_back(const Point2D& point) noexcept {
    push_back(point.x, point.y, point.xerr, point.yerr);
}