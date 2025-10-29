// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <dataset/Dataset2D.h>
#include <utility/Exceptions.h>
#include <dataset/DatasetFactory.h>

#include <vector>

using namespace ausaxs;

Dataset2D::Dataset2D() noexcept : SimpleDataset(0, 4) {}

Dataset2D::Dataset2D(unsigned int rows) noexcept : SimpleDataset(rows, 4) {}

Dataset2D::Dataset2D(std::vector<double> x, std::vector<double> y) noexcept : Dataset2D(x.size()) {
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], 0, 0};
    }
}

Dataset2D::Dataset2D(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) noexcept : Dataset2D(x.size()) {
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], yerr[i], 0};
    }
}

Dataset2D::Dataset2D(std::vector<double> x, std::vector<double> y, std::vector<double> xerr, std::vector<double> yerr) noexcept : Dataset2D(x.size()) {
    for (unsigned int i = 0; i < x.size(); i++) {
        row(i) = {x[i], y[i], yerr[i], xerr[i]};
    }
}

Dataset2D::Dataset2D(const SimpleDataset& data) : Dataset2D(data.size()) {
    for (unsigned int i = 0; i < data.size(); i++) {
        row(i) = {data.x(i), data.y(i), data.yerr(i), 0};
    }
}

Dataset2D::Dataset2D(const io::ExistingFile& path) : Dataset2D() {
    auto dataset = factory::DatasetFactory::construct(path, 4);
    this->data = std::move(dataset->data);
    data.N = dataset->data.N;
}

void Dataset2D::scale_errors(double factor) {
    auto xerr = this->xerr();
    std::transform(xerr.begin(), xerr.end(), xerr.begin(), [&factor] (double val) {return factor*val;});
    SimpleDataset::scale_errors(factor);
}

void Dataset2D::push_back(double x, double y, double xerr, double yerr) {
    data.extend(1);
    row(data.N-1) = {x, y, yerr, xerr};
}

void Dataset2D::push_back(double x, double y) {
    push_back(x, y, 0, 0);
}

void Dataset2D::push_back(const Point2D& point) noexcept {
    push_back(point.x, point.y, point.xerr, point.yerr);
}