// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/Axis.h>

#include <utility/Limit.h>

#include <algorithm>
#include <cmath>
#include <ostream>

using namespace ausaxs;

Axis::Axis() noexcept : bins(0), min(0), max(0) {}

Axis::Axis(const Limit& limits, int bins) noexcept : bins(bins), min(limits.min), max(limits.max) {}

Axis& Axis::operator=(std::initializer_list<double> list) noexcept {
    std::vector<double> d = list;
    bins = static_cast<int>(std::round(d[0])); 
    min = d[1];
    max = d[2];
    return *this;
}

std::string Axis::to_string() const noexcept {
    return "Axis: (" + std::to_string(min) + ", " + std::to_string(max) + ") with " + std::to_string(bins) + " bins";
}

bool Axis::operator==(const Axis& rhs) const noexcept = default;

void Axis::resize(int bins) noexcept {
    auto w = width();
    this->bins = bins;
    this->max = min + bins*w;
}

bool Axis::empty() const noexcept {return bins==0;}

Limit Axis::limits() const noexcept {return {min, max};}

int Axis::get_bin(double value) const noexcept {
    if (bins == 0) [[unlikely]] {return 0;}
    if (value <= min) {return 0;}
    if (value >= max) {return bins;}
    return std::floor((value+1e-6-min)/width()); // +1e-6 to avoid flooring floating point errors, and we will likely never have bins this small anyway
}

double Axis::get_bin_value(int bin) const noexcept {
    if (bins == 0) [[unlikely]] {return 0;}
    return min + bin*width();
}

Axis Axis::sub_axis(double vmin, double vmax) const noexcept {
    int min_bin = get_bin(vmin);
    int max_bin = get_bin(vmax);

    double new_min = get_bin_value(min_bin);
    double new_max = get_bin_value(max_bin);
    return {new_min, new_max, max_bin - min_bin};
}

Axis Axis::sub_axis_covering(double vmin, double vmax) const noexcept {
    if (bins == 0) [[unlikely]] {return *this;}
    int min_bin = std::min(get_bin(vmin), bins-1);
    int max_bin = std::clamp(static_cast<int>(std::ceil((vmax-1e-6-min)/width())), min_bin, bins-1); // -1e-6 to avoid ceiling floating point errors
    return {get_bin_value(min_bin), get_bin_value(max_bin+1), max_bin - min_bin + 1};
}

namespace {
    [[maybe_unused]] std::ostream& operator<<(std::ostream& os, const Axis& axis) noexcept {os << axis.to_string(); return os;}
}