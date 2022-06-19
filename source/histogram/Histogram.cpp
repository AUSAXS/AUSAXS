#include <histogram/Histogram.h>

#include <algorithm>

using namespace hist;

Histogram::Histogram(const Vector<double>& p) noexcept : p(p) {}

Histogram::Histogram(const Vector<double>& p, const Axis& axis) noexcept : p(p), axis(axis) {}

Histogram::Histogram(const Axis& axis) noexcept : p(axis.bins), axis(axis) {}

Histogram& Histogram::operator+=(const Histogram& rhs) noexcept {
    p += rhs.p;
    return *this;
}

Histogram& Histogram::operator-=(const Histogram& rhs) noexcept {
    p -= rhs.p;
    return *this;
}

void Histogram::shorten_axis() {
    int max_bin = 10; // minimum size is 10
    for (int i = axis.bins-1; i >= 10; i--) {
        if (p[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p.resize(max_bin);
    double width = axis.width();
    axis = Axis{max_bin, 0, max_bin*width};
}

void Histogram::generate_axis(unsigned int size) {
    Limit limits = span();
    axis.min = limits.min; axis.max = limits.max;
    axis.bins = size;
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

size_t Histogram::size() const noexcept {return p.size();}

void Histogram::set_plot_options(const plots::PlotOptions& options) {plot_options = options;}

void Histogram::add_plot_options(const std::map<std::string, std::any>& options) {plot_options.set(options);}

void Histogram::add_plot_options(std::string style, std::map<std::string, std::any> options) {plot_options.set(style, options);}

void Histogram::add_plot_options(int color, std::map<std::string, std::any> options) {plot_options.set(color, options);}