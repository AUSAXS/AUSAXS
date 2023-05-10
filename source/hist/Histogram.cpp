#include <hist/Histogram.h>
#include <dataset/SimpleDataset.h>

#include <algorithm>
#include <sstream>

using namespace hist;

Histogram::Histogram(const Vector<double>& p) noexcept : p(p) {}

Histogram::Histogram(const Vector<double>& p, const Axis& axis) noexcept : p(p), axis(axis) {}

Histogram::Histogram(const Axis& axis) noexcept : p(axis.bins), axis(axis) {}

Histogram::~Histogram() = default;

Histogram& Histogram::operator+=(const Histogram& rhs) {
    p += rhs.p;
    return *this;
}

Histogram& Histogram::operator-=(const Histogram& rhs) {
    p -= rhs.p;
    return *this;
}

double& Histogram::operator[](const int i) {
    return p[i];
}

double Histogram::operator[](const int i) const {
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
    axis = Axis{int(max_bin), 0, max_bin*width};
}

void Histogram::generate_axis() {
    Limit limits = span();
    axis.min = limits.min; axis.max = limits.max;
    axis.bins = p.size() == 0 ? 100 : p.size();
}

void Histogram::set_axis(const Axis& axis) noexcept {
    this->axis = axis;
}

Limit Histogram::span() const noexcept {
    if (axis.empty()) {
        auto[min, max] = std::minmax_element(p.begin(), p.end());
        return Limit(*min, *max);
    }
    return axis.limits();
}

Limit Histogram::span_positive() const noexcept {
    if (size() == 0) {
        return Limit(0, 0);
    }

    Limit limits(p[0], p[0]);
    for (double val : p) {
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
