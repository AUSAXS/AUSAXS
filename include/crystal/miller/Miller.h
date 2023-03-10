#pragma once

#include <type_traits>
#include <math.h>
#include <math/Vector3.h>

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

        Vector3<double> normalize() const {
            double length = std::sqrt(h*h + k*k + l*l);
            return Vector3<double>(h/length, k/length, l/length);
        }

        double distance(const Miller& other) const {
            return std::sqrt(std::pow(h - other.h, 2) + std::pow(k - other.k, 2) + std::pow(l - other.l, 2));
        }

        double length() const {return std::sqrt(length2());}
        double length2() const {return h*h + k*k + l*l;}

        bool friedel_equivalent(const Miller& other) {
            return (h == -other.h) && (k == -other.k) && (l == -other.l);
        }

        int h, k, l;
    };
}