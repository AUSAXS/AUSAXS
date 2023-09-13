#include <hist/Histogram.h>
#include <dataset/SimpleDataset.h>

#include <algorithm>
#include <sstream>

using namespace hist;

Histogram::Histogram(const Vector<double>& p) noexcept : p(p) {}

Histogram::Histogram(const Vector<double>& p, const Axis& axis) : p(p), axis(axis) {
    #ifdef DEBUG
        if (p.size() != axis.bins) {throw std::invalid_argument("Histogram: Vector and Axis must have the same number of bins.");}
    #endif
}

Histogram::Histogram(std::vector<double>&& p_tot, const Axis& axis) : p(std::move(p_tot)), axis(axis) {
    #ifdef DEBUG
        if (p.size() != axis.bins) {throw std::invalid_argument("Histogram: Vector and Axis must have the same number of bins.");}
    #endif
}

Histogram::Histogram(const Axis& axis) noexcept : p(axis.bins), axis(axis) {}

Histogram::~Histogram() = default;

Histogram& Histogram::operator+=(const Histogram& rhs) {p += rhs.p; return *this;}
Histogram hist::operator+(const Histogram& lhs, const Histogram& rhs) {
    Histogram result(lhs);
    result += rhs;
    return result;
}

Histogram& Histogram::operator-=(const Histogram& rhs) {p -= rhs.p; return *this;}
Histogram hist::operator-(const Histogram& lhs, const Histogram& rhs) {
    Histogram result(lhs);
    result -= rhs;
    return result;
}

Histogram& Histogram::operator*=(double rhs) {p *= rhs; return *this;}
Histogram hist::operator*(const Histogram& lhs, double rhs) {
    Histogram result(lhs);
    result *= rhs;
    return result;
}

double& Histogram::operator[](int i) {
    return p[i];
}

double Histogram::operator[](int i) const {
    return p[i];
}

void Histogram::resize(unsigned int bins) {
    p.resize(bins);
    axis.resize(bins);
}

void Histogram::shorten_axis(unsigned int min_size) {
    unsigned int max_bin = min_size;
    for (unsigned int i = axis.bins-1; i >= min_size; --i) {
        if (p[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p.resize(max_bin);
    double width = axis.width();
    axis = Axis(0, max_bin*width, int(max_bin));
}

// void Histogram::extend_axis(double qmax) {
//     double width = axis.width();
//     unsigned int bins = int(qmax/width);
//     p.resize(bins);
//     axis = Axis(0, bins*width, bins);
// }

void Histogram::generate_axis() {
    Limit lim = limits();
    axis.min = lim.min; axis.max = lim.max;
    axis.bins = p.size() == 0 ? 100 : p.size();
}

void Histogram::set_axis(const Axis& axis) noexcept {
    this->axis = axis;
}

Limit Histogram::limits() const noexcept {
    if (axis.empty()) {
        auto[min, max] = std::minmax_element(p.begin(), p.end());
        return Limit(*min, *max);
    }
    return axis.limits();
}

Limit Histogram::limits_positive() const noexcept {
    if (size() == 0) {
        return Limit(0, 0);
    }

    Limit limits(p[0], p[0]);
    for (const double val : p) {
        if (0 < val) {
            limits.min = std::min(val, limits.min);
        }
        limits.max = std::max(val, limits.max);
    }
    return limits;
}

std::string Histogram::to_string() const noexcept {
    std::stringstream ss;
    auto ax = axis.as_vector();
    for (unsigned int i = 0; i < size(); i++) {
        ss << ax[i] << " " << p[i] << std::endl;
    }
    return ss.str();
}

unsigned int Histogram::size() const noexcept {return p.size();}

SimpleDataset Histogram::as_dataset() const {
    return SimpleDataset(axis.as_vector(), std::vector<double>(p.begin(), p.end()));
}

bool Histogram::operator==(const Histogram& rhs) const {
    return p == rhs.p && axis == rhs.axis;
}