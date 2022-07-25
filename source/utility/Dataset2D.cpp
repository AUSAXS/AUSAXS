#include <utility/Dataset2D.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <utility/Settings.h>

#include <vector>
#include <string>
#include <fstream>
#include <random>

using std::vector, std::string;

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
    push_back(point.x, point.y, point.xerr, point.yerr);
}

void Dataset2D::load(std::string path) {
    Dataset::load(path);
    if (M != 4) {
        throw except::io_error("Error in Dataset2D::load: Dataset has wrong number of columns.");
    }
    limit_x(setting::fit::q_low, setting::fit::q_high);
}