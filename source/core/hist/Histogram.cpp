/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/Histogram.h>
#include <dataset/SimpleDataset.h>

#include <algorithm>
#include <sstream>
#include <cassert>

using namespace hist;

Histogram::Histogram(const Vector<double>& p) noexcept : p(p) {generate_axis();}
Histogram::Histogram(const Vector<double>& p, const Axis& axis) : p(p), axis(axis) {}
Histogram::Histogram(std::vector<double>&& p_tot, const Axis& axis) : p(std::move(p_tot)), axis(axis) {}
Histogram::Histogram(const Axis& axis) noexcept : p(axis.bins), axis(axis) {}

Histogram& Histogram::operator+=(const Histogram& rhs) {p += rhs.p; return *this;}
Histogram hist::operator+(const Histogram& lhs, const Histogram& rhs) {
    Histogram result(lhs);
    result += rhs;
    return result;
}

Histogram& Histogram::operator-=(const Histogram& rhs) {p -= rhs.p; return *this;}
Histogram hist::operator-(const Histogram& lhs, const Histogram& rhs) {
    Histogram result(lhs);
    result -= rhs;
    return result;
}

Histogram& Histogram::operator*=(double rhs) {p *= rhs; return *this;}
Histogram hist::operator*(const Histogram& lhs, double rhs) {
    Histogram result(lhs);
    result *= rhs;
    return result;
}

double& Histogram::operator[](int i) {
    return p[i];
}

double Histogram::operator[](int i) const {
    return p[i];
}

void Histogram::resize(unsigned int bins) {
    p.resize(bins);
    axis.resize(bins);
}

void Histogram::shorten_axis(unsigned int min_size) {
    if (p.size() < min_size) {return;}
    unsigned int max_bin = min_size;
    for (unsigned int i = p.size()-1; i > min_size; --i) {
        if (p[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p.resize(max_bin);
    if (!axis.empty()) {
        axis = Axis(axis.min, axis.min + max_bin*axis.width(), max_bin);
    }
}

const Axis& Histogram::get_axis() const {
    return axis;
}

const std::vector<double>& Histogram::get_counts() const {
    return p.data;
}

std::vector<double>& Histogram::get_counts() {
    return const_cast<std::vector<double>&>(const_cast<const Histogram*>(this)->get_counts());
}

double& Histogram::get_count(unsigned int i) {
    assert(i < p.size() && "Index out of bounds");
    return p[i];
}

const double& Histogram::get_count(unsigned int i) const {
    assert(i < p.size() && "Index out of bounds");
    return p[i];
}

double& Histogram::index(unsigned int i) {
    return get_count(i);
}

const double& Histogram::index(unsigned int i) const {
    return get_count(i);
}

void Histogram::generate_axis() {
    axis = Axis(0, p.size(), p.size());
}

void Histogram::set_axis(const Axis& axis) noexcept {
    this->axis = axis;
}

void Histogram::normalize(double sum) {
    double total = std::accumulate(p.begin(), p.end(), 0.0);
    assert(total != 0 && "Cannot normalize a histogram with a sum of 0");
    std::transform(p.begin(), p.end(), p.begin(), [sum, total] (double x) {return x/total*sum;});
}

void Histogram::normalize_max(double max) {
    double max_val = *std::max_element(p.begin(), p.end());
    assert(max_val != 0 && "Cannot normalize a histogram with a max of 0");
    std::transform(p.begin(), p.end(), p.begin(), [max_val, max] (double x) {return x/max_val*max;});
}

void Histogram::add_count(unsigned int i, double count) {
    assert(i < p.size() && "Index out of bounds");
    p[i] += count;
}

void Histogram::set_count(unsigned int i, double count) {
    assert(i < p.size() && "Index out of bounds");
    p[i] = count;
}

void Histogram::set_count(const std::vector<double>& counts) {
    p = counts;
}

void Histogram::bin(double value) {
    p[std::min(axis.get_bin(value), p.size()-1)]++;
}

void Histogram::bin(const std::vector<double>& values) {
    std::for_each(values.begin(), values.end(), [this] (double v) {bin(v);});
}

Limit Histogram::span_y() const noexcept {
    if (p.size() == 0) {return Limit(0, 0);}
    auto[min, max] = std::minmax_element(p.begin(), p.end());
    return Limit(*min, *max);
}

Limit Histogram::span_y_positive() const noexcept {
    if (size() == 0) {
        return Limit(0, 0);
    }

    double min = 0;
    double max = 0;
    for (unsigned int i = 0; i < p.size(); ++i) {
        if (0 < p[i]) {
            min = std::min(p[i], min);
        }
        max = std::max(p[i], max);
    }
    return Limit(min, max);
}

std::string Histogram::to_string() const noexcept {
    std::stringstream ss;
    auto ax = axis.as_vector();
    for (unsigned int i = 0; i < size(); i++) {
        ss << ax[i] << " " << p[i] << std::endl;
    }
    return ss.str();
}

unsigned int Histogram::size() const noexcept {return p.size();}

SimpleDataset Histogram::as_dataset() const {
    return SimpleDataset(axis.as_vector(), std::vector<double>(p.begin(), p.end()));
}

bool Histogram::operator==(const Histogram& rhs) const {
    return p == rhs.p && axis == rhs.axis;
}