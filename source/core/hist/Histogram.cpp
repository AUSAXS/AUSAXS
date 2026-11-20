// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/Histogram.h>

#include <dataset/Dataset.h>

#include <algorithm>
#include <cassert>
#include <numeric>
#include <sstream>

using namespace ausaxs;
using namespace ausaxs::hist;

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

void Histogram::resize(int bins) {
    p.resize(bins);
    axis.resize(bins);
}

const Axis& Histogram::get_axis() const {
    return axis;
}

const std::vector<double>& Histogram::get_counts() const {
    return p.data;
}

double& Histogram::get_count(int i) {
    assert(i < p.size() && "Index out of bounds");
    return p[i];
}

const double& Histogram::get_count(int i) const {
    assert(i < p.size() && "Index out of bounds");
    return p[i];
}

double& Histogram::index(int i) {
    return get_count(i);
}

const double& Histogram::index(int i) const {
    return get_count(i);
}

void Histogram::generate_axis() {
    axis = Axis(0, p.size(), p.size());
}

void Histogram::set_axis(const Axis& axis) noexcept {
    this->axis = axis;
}

void Histogram::normalize(double sum) {
    double total = std::reduce(p.begin(), p.end(), 0.0);
    assert(total != 0 && "Cannot normalize a histogram with a sum of 0");
    std::ranges::transform(p, p.begin(), [sum, total] (double x) {return x/total*sum;});
}

void Histogram::normalize_max(double max) {
    double max_val = *std::ranges::max_element(p);
    assert(max_val != 0 && "Cannot normalize a histogram with a max of 0");
    std::ranges::transform(p, p.begin(), [max_val, max] (double x) {return x/max_val*max;});
}

void Histogram::add_count(int i, double count) {
    assert(i < p.size() && "Index out of bounds");
    p[i] += count;
}

void Histogram::set_count(int i, double count) {
    assert(i < p.size() && "Index out of bounds");
    p[i] = count;
}

void Histogram::set_count(const std::vector<double>& counts) {
    p = counts;
}

Limit Histogram::span_y() const noexcept {
    if (p.empty()) {return {0, 0};}
    auto[min, max] = std::ranges::minmax_element(p);
    return {*min, *max};
}

std::string Histogram::to_string() const noexcept {
    std::stringstream ss;
    auto ax = axis.as_vector();
    for (int i = 0; i < size(); i++) {
        ss << ax[i] << " " << p[i] << std::endl;
    }
    return ss.str();
}

int Histogram::size() const noexcept {return p.size();}

Dataset Histogram::as_dataset() const {
    return {axis.as_vector(), p.data};
}

bool Histogram::operator==(const Histogram& rhs) const {
    return p == rhs.p && axis == rhs.axis;
}