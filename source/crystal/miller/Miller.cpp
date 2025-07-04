// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <crystal/miller/Miller.h>

#include <math/Vector3.h>

using namespace ausaxs::crystal;

Miller::Miller() : h(0), k(0), l(0) {}
Miller::Miller(int h, int k, int l) : h(h), k(k), l(l) {}

bool Miller::operator==(const Miller& other) const {return h == other.h && k == other.k && l == other.l;}
bool Miller::operator!=(const Miller& other) const {return !(*this == other);}
Miller Miller::operator*(int n) const {return Miller(h*n, k*n, l*n);}

ausaxs::Vector3<double> Miller::normalize() const {
    double length = std::sqrt(h*h + k*k + l*l);
    return Vector3<double>(h/length, k/length, l/length);
}

double Miller::distance(const Miller& other) const {
    return std::sqrt(std::pow(h - other.h, 2) + std::pow(k - other.k, 2) + std::pow(l - other.l, 2));
}

double Miller::length() const {return std::sqrt(length2());}
double Miller::length2() const {return h*h + k*k + l*l;}

bool Miller::friedel_equivalent(const Miller& other) {
    return (h == -other.h) && (k == -other.k) && (l == -other.l);
}