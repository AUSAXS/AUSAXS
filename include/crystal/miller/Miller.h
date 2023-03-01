#pragma once

#include <type_traits>
#include <math.h>

namespace crystal {
    /**
     * @brief A struct to represent a Miller index. 
     */
    struct Miller {
        Miller() : h(0), k(0), l(0) {}
        Miller(int h, int k, int l) : h(h), k(k), l(l) {}

        bool operator==(const Miller& other) const {return h == other.h && k == other.k && l == other.l;}
        bool operator!=(const Miller& other) const {return !(*this == other);}
        Miller operator*(int n) const {return Miller(h*n, k*n, l*n);}

        double length() const {return std::sqrt(length2());}
        double length2() const {return h*h + k*k + l*l;}

        bool friedel_equivalent(const Miller& other) {
            return (h == -other.h) && (k == -other.k) && (l == -other.l);
        }

        int h, k, l;
    };
}