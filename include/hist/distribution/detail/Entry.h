#pragma once

#include <constants/Axes.h>

namespace hist {
    namespace detail {
        struct Entry {
            Entry() = default;
            Entry(constants::axes::d_type count) : count(count) {}
            Entry(constants::axes::d_type count, double content) : count(count), content(content) {}

            double count = 0;
            double content = 0;

            /**
             * @brief Add the distance to this bin, and increase the counter by one.
             */
            void add(double distance) {
                ++count;
                content += distance;
            }

            int operator+(const Entry& other) const {
                return count + other.count;
            }

            Entry operator+(const Entry& other) {
                return Entry(count + other.count, content + other.content);
            }

            Entry& operator+=(const Entry& other) {
                count += other.count;
                content += other.content;
                return *this;
            }

            int operator-(const Entry& other) const {
                return count - other.count;
            }

            Entry operator-(const Entry& other) {
                return Entry(count - other.count, content - other.content);
            }

            Entry& operator-=(const Entry& other) {
                count -= other.count;
                content -= other.content;
                return *this;
            }

            friend Entry operator*(const Entry& entry, double factor) {
                return Entry(entry.count*factor, entry.content*factor);
            }

            friend Entry operator*(double factor, const Entry& entry) {
                return Entry(entry.count*factor, entry.content*factor);
            }

            bool operator==(double other) const {
                return count == other;
            }
        };
    }
}