#include <utility/Dataset2D.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <utility/Settings.h>

#include <vector>
#include <string>
#include <fstream>
#include <random>

using std::vector, std::string;

Dataset2D::Dataset2D() noexcept : SimpleDataset(0, 4) {}

Dataset2D::Dataset2D(unsigned int rows) noexcept : SimpleDataset(rows, 4) {}

Dataset2D::Dataset2D(std::vector<double> x, std::vector<double> y) noexcept : Dataset2D(x.size()) {
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], 0, 0};
    }
}

Dataset2D::Dataset2D(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel) : Dataset2D(x, y) {
    set_col_names({xlabel, ylabel, std::string(ylabel)+"err", std::string(xlabel)+"err"});
    options.xlabel = xlabel;
    options.ylabel = ylabel;
}

Dataset2D::Dataset2D(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) noexcept : Dataset2D(x.size()) {
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], yerr[i], 0};
    }
}

Dataset2D::Dataset2D(std::vector<double> x, std::vector<double> y, std::vector<double> xerr, std::vector<double> yerr) noexcept : Dataset2D(x.size()) {
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], xerr[i], yerr[i]};
    }
}

Dataset2D::Dataset2D(const SimpleDataset& data) : Dataset2D(data.size()) {
    for (unsigned int i = 0; i < data.size(); i++) {
        row(i) = {data.x(i), data.y(i), data.yerr(i), 0};
    }
}

Dataset2D::Dataset2D(std::string path) : Dataset2D() {
    load(path);
}

void Dataset2D::scale_errors(double factor) {
    auto xerr = this->xerr();
    auto yerr = this->yerr();
    std::transform(xerr.begin(), xerr.end(), xerr.begin(), [&factor] (double val) {return factor*val;});
    std::transform(yerr.begin(), yerr.end(), yerr.begin(), [&factor] (double val) {return factor*val;});
}

void Dataset2D::push_back(double x, double y, double xerr, double yerr) {
    extend(1);
    row(N-1) = {x, y, yerr, xerr};
}

void Dataset2D::push_back(double x, double y) {
    push_back(x, y, 0, 0);
}

void Dataset2D::push_back(const Point2D& point) noexcept {
    push_back(point.x, point.y, point.yerr, point.xerr);
}

void Dataset2D::load(std::string path) {
    Dataset::load(path);
    if (M != 4) {
        throw except::io_error("Error in Dataset2D::load: Dataset has wrong number of columns.");
    }
    names = {"q", "I", "Ierr", "qerr"}; // set column names
    unsigned int N = size();
    limit_x(setting::axes::qmin, setting::axes::qmax);
    if (N != size() && setting::general::verbose) {
        std::cout << "\tRemoved " << N - size() << " data points outside specified q-range [" << setting::axes::qmin << ", " << setting::axes::qmax << "]." << std::endl;
    }
}