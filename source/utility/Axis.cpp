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

std::ostream& operator<<(std::ostream& os, const Axis& axis) noexcept {os << axis.to_string(); return os;}