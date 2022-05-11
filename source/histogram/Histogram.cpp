#include <histogram/Histogram.h>

#include <algorithm>

using namespace hist;

Histogram::Histogram(const Vector<double>& p) : p(p) {}

Histogram::Histogram(const Vector<double>& p, const Axis& axis) : p(p), axis(axis) {}

Histogram::Histogram(const Axis& axis) : p(axis.bins), axis(axis) {}

Histogram& Histogram::operator+=(const Histogram& rhs) {
    p += rhs.p;
    return *this;
}

Histogram& Histogram::operator-=(const Histogram& rhs) {
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
    auto[min, max] = std::minmax_element(p.begin(), p.end());
    axis.min = *min; axis.max = *max;
    axis.bins = size;
}

void Histogram::set_axis(const Axis& axis) {
    this->axis = axis;
}

size_t Histogram::size() const {return p.size();}

void Histogram::set_plot_options(const plots::PlotOptions& options) {plot_options = options;}

void Histogram::add_plot_options(const std::map<std::string, std::any>& options) {plot_options.set(options);}

void Histogram::add_plot_options(std::string style, std::map<std::string, std::any> options) {plot_options.set(style, options);}

void Histogram::add_plot_options(int color, std::map<std::string, std::any> options) {plot_options.set(color, options);}