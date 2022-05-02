#include <histogram/Histogram.h>

#include <algorithm>

using namespace histogram;

Histogram::Histogram(const Vector<double>& p) : p(p) {}

Histogram::Histogram(const Vector<double>& p, const Axis& axis) : p(p), axis(axis) {}

Histogram::Histogram(const Axis& axis) : p(axis.bins), axis(axis) {}

Histogram& Histogram::operator+=(const Histogram& rhs) {
    p += rhs.p;
    return *this;
}

Histogram& Histogram::operator-=(const Histogram& rhs) {
    p -= rhs.p;
    return *this;
}

void Histogram::shorten_axis() {
    int max_bin = 10; // minimum size is 10
    for (int i = axis.bins-1; i >= 10; i--) {
        if (p[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p.resize(max_bin);
    double width = axis.width();
    axis = Axis{max_bin, 0, max_bin*width};
}