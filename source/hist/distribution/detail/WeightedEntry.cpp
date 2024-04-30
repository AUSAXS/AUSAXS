/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/distribution/detail/WeightedEntry.h>

#include <iostream>

using namespace hist::detail;

WeightedEntry::WeightedEntry() = default;
WeightedEntry::WeightedEntry(constants::axes::d_type value, unsigned int count, double bin_center) : value(value), count(count), bin_center(bin_center) {}

void WeightedEntry::add(float distance, double value) {
    ++count;
    bin_center += distance;
    this->value += value;
}

void WeightedEntry::add2(float distance, double value) {
    count += 2;
    bin_center += 2*distance;
    this->value += 2*value;
}

WeightedEntry WeightedEntry::operator+(const WeightedEntry& other) const {
    return WeightedEntry(value + other.value, count + other.count, bin_center + other.bin_center);
}

WeightedEntry& WeightedEntry::operator+=(const WeightedEntry& other) {
    count += other.count;
    value += other.value;
    bin_center += other.bin_center;
    return *this;
}

WeightedEntry WeightedEntry::operator-(const WeightedEntry& other) const {
    return WeightedEntry(value - other.value, count - other.count, bin_center - other.bin_center);
}

WeightedEntry& WeightedEntry::operator-=(const WeightedEntry& other) {
    count -= other.count;
    value -= other.value;
    bin_center -= other.bin_center;
    return *this;
}

bool WeightedEntry::operator==(double other) const {
    return value == other;
}

WeightedEntry hist::detail::operator*(const WeightedEntry& entry, double factor) {
    return WeightedEntry(entry.value*factor, factor*entry.count, entry.bin_center*factor);
}

WeightedEntry hist::detail::operator*(double factor, const WeightedEntry& entry) {
    return entry*factor;
}

std::ostream& hist::detail::operator<<(std::ostream& os, const WeightedEntry& entry) {
    os << "WeightedEntry(" << entry.value << ", " << entry.count << ", " << entry.bin_center << ")";
    return os;
}