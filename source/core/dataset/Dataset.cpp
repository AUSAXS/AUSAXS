// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <dataset/Dataset.h>
#include <math/MovingAverager.h>
#include <math/CubicSpline.h>
#include <math/PeakFinder.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <utility/Console.h>
#include <dataset/DatasetFactory.h>
#include <settings/GeneralSettings.h>

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using namespace ausaxs;

Dataset::Dataset() = default;
Dataset::Dataset(const Dataset& d) = default;
Dataset::Dataset(Dataset&& d) = default;
Dataset& Dataset::operator=(const Dataset& other) = default;
Dataset& Dataset::operator=(Dataset&& other) = default;

Dataset::Dataset(Matrix<double>&& m) : data(std::move(m)) {
    set_default_names();
}

Dataset::Dataset(const std::vector<std::string>& col_names) : data(0, col_names.size()), names(col_names) {}

Dataset::Dataset(const std::vector<std::vector<double>>& cols, const std::vector<std::string>& col_names) : data(cols), names(col_names) {}

Dataset::Dataset(unsigned int rows, unsigned int cols) : data(rows, cols) {
    set_default_names();
}

Dataset::Dataset(const std::vector<std::vector<double>>& cols) : data(cols) {
    set_default_names();
}

Dataset::Dataset(const io::ExistingFile& path) : Dataset() {
    *this = std::move(*factory::DatasetFactory::construct(path));
    set_default_names();
}

Dataset::~Dataset() = default;

void Dataset::assign_matrix(Matrix<double>&& m) {
    if (m.M != data.M) {
        throw except::invalid_operation("Dataset::operator=: Matrix has wrong number of columns. "
        "Expected " + std::to_string(data.M) + ", but got " + std::to_string(m.M));
    }
    force_assign_matrix(std::move(m));
}

void Dataset::force_assign_matrix(Matrix<double>&& m) {
    data.data = std::move(m.data);
    data.N = m.N;
    data.M = m.M;
}

bool Dataset::empty() const noexcept {
    return data.data.empty();
}

void Dataset::limit_x(const Limit& limits) {
    if (size() == 0) {return;}
    if (limits.min < x(0) && x(size()-1) < limits.max) {return;}

    Matrix<double> limited(0, data.M); 
    for (unsigned int i = 0; i < size(); i++) {
        double val = x(i);
        if (val < limits.min) {continue;}
        else if (limits.max < val) {break;}
        limited.push_back(row(i));
    }
    assign_matrix(std::move(limited));
}

void Dataset::limit_y(const Limit& limits) {
    if (size() == 0) {return;}

    Matrix<double> limited(0, data.M);
    for (unsigned int i = 0; i < size(); i++) {
        double val = y(i);
        if (val < limits.min || limits.max < val) {continue;}
        limited.push_back(row(i));
    }
    assign_matrix(std::move(limited));
}

void Dataset::limit_x(double min, double max) {limit_x({min, max});}
void Dataset::limit_y(double min, double max) {limit_y({min, max});}

MutableColumn<double> Dataset::col(std::string_view column) {
    for (unsigned int i = 0; i < names.size(); ++i) {
        if (names[i] == column) {
            return col(i);
        }
    }
    throw except::invalid_operation("Dataset::col: Column \"" + std::string(column) + "\" not found. Available columns:\n\t" + utility::join(names, "\n\t"));
}

const ConstColumn<double> Dataset::col(std::string_view column) const {
    for (unsigned int i = 0; i < names.size(); ++i) {
        if (names[i] == column) {
            return col(i);
        }
    }
    throw except::invalid_operation("Dataset::col: Column \"" + std::string(column) + "\" not found. Available columns: " + utility::join(names, "\n"));
}

MutableColumn<double> Dataset::col(unsigned int index) {
    return data.col(index);
}

const ConstColumn<double> Dataset::col(unsigned int index) const {
    return data.col(index);
}

MutableRow<double> Dataset::row(unsigned int index) {
    return data.row(index);
}

const ConstRow<double> Dataset::row(unsigned int index) const {
    return data.row(index);
}

void Dataset::set_col_names(const std::vector<std::string>& names) {
    if (names.size() != data.M) {
        throw except::invalid_operation(
            "Dataset::set_col_names: Number of names does not match number of columns. "
            "(" + std::to_string(names.size()) + " != " + std::to_string(data.M) + ")"
        );
    }
    this->names = names;
}

void Dataset::set_col_names(unsigned int i, const std::string& name) {
    names[i] = name;
}

std::vector<std::string> Dataset::get_col_names() {
    return names;
}

std::string Dataset::get_col_names(unsigned int i) {
    return names[i];
}

Dataset Dataset::select_columns(const std::vector<unsigned int>& cols) const {
    if (cols.size() == 0) {throw except::invalid_argument("Dataset::select_columns: No columns selected.");}
    Dataset data(this->data.N, cols.size());
    std::vector<std::string> col_names(cols.size());
    for (unsigned int i = 0; i < cols.size(); i++) {
        data.col(i) = col(cols[i]);
        col_names[i] = names[cols[i]];
    }
    data.set_col_names(col_names);
    return data;
}

Dataset Dataset::select_columns(const std::vector<std::string>& cols) const {
    std::vector<unsigned int> indices(cols.size());
    std::transform(cols.begin(), cols.end(), indices.begin(), [this](const std::string& name) {
        for (unsigned int i = 0; i < names.size(); i++) {
            if (names[i] == name) {return i;}
        }
        throw except::invalid_argument("Dataset::select_columns: Column \"" + name + "\" not found.");
    });
    return select_columns(indices);
}

void Dataset::save(const io::File& path, const std::string& header) const {
    path.directory().create();

    // check if file was succesfully opened
    std::ofstream output(path);
    if (!output.is_open()) {throw std::ios_base::failure("IntensityFitter::save: Could not open file \"" + path.str() + "\"");}

    // write header
    if (!header.empty()) {
        output << header << std::endl;
    }

    // write column titles
    if (names.size() < data.M) {
        throw except::unexpected(
            "Dataset::save: Number of column names (" + std::to_string(names.size()) + ") "
            "does not match number of columns (" + std::to_string(data.M) + ")."
        );
    }
    for (unsigned int j = 0; j < data.M; j++) {
        output << std::left << std::setw(16) << names[j] << "\t";
    }
    output << std::endl;

    // write data
    for (unsigned int i = 0; i < data.N; i++) {
        for (unsigned int j = 0; j < data.M-1; j++) {
            output << std::left << std::setw(16) << std::setprecision(8) << std::scientific << index(i, j) << "\t";
        }
        output << index(i, data.M-1) << "\n";
    }
    output.close();
}

void Dataset::load(const io::ExistingFile& path) {
    auto dataset = factory::DatasetFactory::construct(path, data.M);
    if (dataset->data.M != data.M) {
        throw except::invalid_operation(
            "Dataset::load: Number of columns does not match. "
            "(" + std::to_string(dataset->data.M) + " != " + std::to_string(data.M) + ")"
        );
    }
    *this = std::move(*dataset);
    set_default_names();
}

void Dataset::set_default_names() {
    names.resize(data.M);
    for (unsigned int i = 0; i < data.M; i++) {
        names[i] = "col_" + std::to_string(i);
    }
}

Dataset Dataset::rolling_average(unsigned int window_size) const {
    Dataset result(*this);
    result.y() = MovingAverage::average_half(y(), window_size);
    return result;
}

Dataset Dataset::interpolate(unsigned int n) const {
    Matrix<double> interpolated(size()*(n+1)-n-1, data.M);

    std::vector<math::CubicSpline> splines;
    for (unsigned int col_index = 1; col_index < data.M; ++col_index) {
        splines.push_back(math::CubicSpline(x(), col(col_index)));
    }

    for (unsigned int i = 0; i < size()-1; i++) {
        double x = this->x(i);
        interpolated[i*(n+1)] = row(i);

        double x_next = this->x(i+1);
        double step = (x_next - x)/(n+1);
        for (unsigned int j = 0; j < n; j++) {
            std::vector<double> row_new(data.M);
            row_new[0] = x + (j+1)*step;;
            for (unsigned int k = 1; k < data.M; k++) {
                row_new[k] = splines[k-1].spline(row_new[0]);
            }
            interpolated[i*(n+1) + j + 1] = row_new;
        }
    }
    return interpolated;
}

std::vector<double> Dataset::find_minimum(unsigned int col_i) const {
    if (size() == 0) {
        if (settings::general::verbose) {
            console::print_warning("Warning in Dataset::find_minimum: Dataset is empty.");
        }
        return std::vector<double>(data.M, 0);
    }
    
    unsigned int min_index = 0;
    double min_value = y(0);
    for (unsigned int i = 1; i < size(); i++) {
        if (col(col_i)[i] < min_value) {
            min_index = i;
            min_value = col(col_i)[i];
        }
    }
    return row(min_index);
}

Dataset Dataset::interpolate(const std::vector<double>& newx) const {
    Matrix<double> interpolated(newx.size(), data.M);

    std::vector<math::CubicSpline> splines;
    for (unsigned int col_index = 1; col_index < data.M; ++col_index) {
        splines.push_back(math::CubicSpline(x(), col(col_index)));
    }

    for (unsigned int i = 0; i < newx.size(); i++) {
        std::vector<double> row_new(data.M);
        row_new[0] = newx[i];
        for (unsigned int j = 0; j < splines.size(); j++) {
            row_new[1+j] = splines[j].spline(newx[i]);
        }
        interpolated[i] = row_new;
    }
    return interpolated;
}

double Dataset::interpolate_x(double x, unsigned int col_index) const {
    math::CubicSpline spline(this->x(), col(col_index));
    return spline.spline(x);
}

void Dataset::append(const Dataset& other) {
    if (this == &other) {throw except::invalid_argument("Dataset::append: Cannot append to itself.");}
    if (data.M != other.data.M) {throw except::invalid_argument("Dataset::append: Number of columns does not match.");}
    unsigned int n = size();
    data.extend(other.size());
    for (unsigned int i = 0; i < other.size(); i++) {
        row(n+i) = other.row(i);
    }
}

void Dataset::sort_x() {
    Matrix<double> newdata(data.N, data.M);
    std::vector<unsigned int> indices(data.N);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [this] (unsigned int i, unsigned int j) {return x(i) < x(j);});
    for (unsigned int i = 0; i < data.N; i++) {
        newdata.row(i) = this->row(indices[i]);
    }
    this->assign_matrix(std::move(newdata));
}

std::string Dataset::to_string() const {
    std::stringstream ss;
    for (unsigned int i = 0; i < size(); i++) {
        for (unsigned int j = 0; j < data.M; j++) {
            ss << std::setw(16) << std::setprecision(8) << std::scientific << index(i, j) << " ";
        }
        ss << "\n";
    }
    return ss.str();
}

double Dataset::index(unsigned int i, unsigned int j) const {
    return data.index(i, j);
}

double& Dataset::index(unsigned int i, unsigned int j) {
    return data.index(i, j);
}

void Dataset::push_back(const std::vector<double>& row) {
    data.push_back(row);
}

std::vector<unsigned int> Dataset::find_minima(unsigned int min_spacing, double min_prominence) const {
    return math::find_minima(x(), y(), min_spacing, min_prominence);
}

std::vector<unsigned int> Dataset::find_maxima(unsigned int min_spacing, double min_prominence) const {
    return math::find_minima(x(), -y(), min_spacing, min_prominence);
}

bool Dataset::is_named() const noexcept {
    for (unsigned int i = 0; i < data.M; i++) {
        if (names[i] != "col_" + std::to_string(i)) {
            return true;
        }
    }
    return false;
}

unsigned int Dataset::size() const noexcept {
    return size_rows();
}

unsigned int Dataset::size_rows() const noexcept {
    return data.N;
}

unsigned int Dataset::size_cols() const noexcept {
    return data.M;
}

bool Dataset::operator==(const Dataset& other) const = default;