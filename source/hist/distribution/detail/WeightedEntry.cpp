#include <hist/distribution/detail/WeightedEntry.h>

using namespace hist::detail;

WeightedEntry::WeightedEntry() = default;
WeightedEntry::WeightedEntry(constants::axes::d_type value, unsigned int count, double bin_center) : value(value), count(count), bin_center(bin_center) {}

void WeightedEntry::add(float distance, double value) {
    ++count;
    bin_center += distance;
    this->value += value;
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