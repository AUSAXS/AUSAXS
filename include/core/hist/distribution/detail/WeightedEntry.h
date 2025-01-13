#pragma once

#include <constants/ConstantsAxes.h>
#include <utility/TypeTraits.h>

#include <iosfwd>

namespace ausaxs::hist::detail {
    /**
     * @brief This struct is a small wrapper around a binned value, automatically keeping tracking the weighted center of the bin. 
     *        This is used to keep the value and weight tracking variables close together in memory, thus improving cache locality for the distance calculations.
     */
    struct WeightedEntry {
        WeightedEntry();
        WeightedEntry(constants::axes::d_type value, unsigned int count, double bin_center);

        /**
         * @brief Add the distance to this bin, and increase the counter by one.
         */
        template<int N>
        void add(float distance, double value);

        WeightedEntry operator+(const WeightedEntry& other) const;
        WeightedEntry& operator+=(const WeightedEntry& other);
        WeightedEntry operator-(const WeightedEntry& other) const;
        WeightedEntry& operator-=(const WeightedEntry& other);
        bool operator==(double other) const;
        
        constants::axes::d_type value = 0;
        unsigned int count = 0;
        double bin_center = 0;
    };

    WeightedEntry operator*(const WeightedEntry& entry, double factor);
    WeightedEntry operator*(double factor, const WeightedEntry& entry);
    std::ostream& operator<<(std::ostream& os, const WeightedEntry& entry);
}
static_assert(supports_nothrow_move_v<ausaxs::hist::detail::WeightedEntry>, "WeightedEntry should support nothrow move semantics.");



inline ausaxs::hist::detail::WeightedEntry::WeightedEntry() = default;
inline ausaxs::hist::detail::WeightedEntry::WeightedEntry(constants::axes::d_type value, unsigned int count, double bin_center) : value(value), count(count), bin_center(bin_center) {}

template<int N>
inline void ausaxs::hist::detail::WeightedEntry::add(float distance, double value) {
    count += N;
    bin_center += N*distance;
    this->value += N*value;
}

inline ausaxs::hist::detail::WeightedEntry ausaxs::hist::detail::WeightedEntry::operator+(const WeightedEntry& other) const {
    return WeightedEntry(value + other.value, count + other.count, bin_center + other.bin_center);
}

inline ausaxs::hist::detail::WeightedEntry& ausaxs::hist::detail::WeightedEntry::operator+=(const WeightedEntry& other) {
    count += other.count;
    value += other.value;
    bin_center += other.bin_center;
    return *this;
}

inline ausaxs::hist::detail::WeightedEntry ausaxs::hist::detail::WeightedEntry::operator-(const WeightedEntry& other) const {
    return WeightedEntry(value - other.value, count - other.count, bin_center - other.bin_center);
}

inline ausaxs::hist::detail::WeightedEntry& ausaxs::hist::detail::WeightedEntry::operator-=(const WeightedEntry& other) {
    count -= other.count;
    value -= other.value;
    bin_center -= other.bin_center;
    return *this;
}

inline bool ausaxs::hist::detail::WeightedEntry::operator==(double other) const {
    return value == other;
}

inline ausaxs::hist::detail::WeightedEntry ausaxs::hist::detail::operator*(const WeightedEntry& entry, double factor) {
    return WeightedEntry(entry.value*factor, factor*entry.count, entry.bin_center*factor);
}

inline ausaxs::hist::detail::WeightedEntry ausaxs::hist::detail::operator*(double factor, const WeightedEntry& entry) {
    return entry*factor;
}