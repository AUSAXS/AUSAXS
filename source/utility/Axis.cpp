#include <utility/Axis.h>
#include <utility/Limit.h>

#include <ostream>

Axis::Axis() noexcept : bins(0), min(0), max(0) {}

Axis::Axis(const Limit& limits, int bins) noexcept : bins(bins), min(limits.min), max(limits.max) {}

Axis& Axis::operator=(std::initializer_list<double> list) noexcept {
    std::vector<double> d = list;
    bins = std::round(d[0]); 
    min = d[1];
    max = d[2];
    return *this;
}

std::string Axis::to_string() const noexcept {
    return "Axis: (" + std::to_string(min) + ", " + std::to_string(max) + ") with " + std::to_string(bins) + " bins";
}

bool Axis::operator==(const Axis& rhs) const noexcept {
    if (bins != rhs.bins) {return false;}
    if (min != rhs.min) {return false;}
    if (max != rhs.max) {return false;}
    return true;
}

bool Axis::operator!=(const Axis& rhs) const noexcept {return !operator==(rhs);}

void Axis::resize(unsigned int bins) noexcept {
    auto w = width();
    this->bins = bins;
    this->max = min + bins*w;
}

bool Axis::empty() const noexcept {return bins==0;}

Limit Axis::limits() const noexcept {return Limit(min, max);}

unsigned int Axis::get_bin(double value) const noexcept {
    if (bins == 0) [[unlikely]] {return 0;}
    if (value <= min) {return 0;}
    if (value >= max) {return bins;}
    return std::floor((value+1e-6-min)/width()); // +1e-6 to avoid flooring floating point errors, and we will likely never have bins this small anyway
}

Axis Axis::sub_axis(double vmin, double vmax) const noexcept {
    unsigned int min_bin = get_bin(vmin);
    unsigned int max_bin = get_bin(vmax);
    double width = this->width();

    double new_min = this->min + min_bin*width;
    double new_max = this->min + max_bin*width;

    return Axis(new_min, new_max, max_bin - min_bin);
}

std::ostream& operator<<(std::ostream& os, const Axis& axis) noexcept {os << axis.to_string(); return os;}