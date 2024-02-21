#include <math/SimpleLeastSquares.h>
#include <math/Statistics.h>
#include <math/MovingAverager.h>
#include <math/CubicSpline.h>
#include <math/PeakFinder.h>
#include <dataset/Dataset.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <utility/Console.h>
#include <dataset/DatasetFactory.h>
#include <dataset/SimpleDataset.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>
#include <settings/GeneralSettings.h>

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

Dataset::Dataset() = default;

Dataset::Dataset(const Dataset& d) = default;

Dataset::Dataset(Dataset&& d) = default;

Dataset::Dataset(Matrix&& m) : Matrix(std::move(m)) {
    set_default_names();
}

Dataset::Dataset(const std::vector<std::string>& col_names) : Matrix(0, col_names.size()), names(col_names) {}

Dataset::Dataset(const std::vector<std::vector<double>>& cols, const std::vector<std::string>& col_names) : Matrix(cols), names(col_names) {}

Dataset::Dataset(unsigned int rows, unsigned int cols) : Matrix(rows, cols) {
    set_default_names();
}

Dataset::Dataset(const std::vector<std::vector<double>>& cols) : Matrix(cols) {
    set_default_names();
}

Dataset::Dataset(const io::ExistingFile& path) : Dataset() {
    *this = std::move(*factory::DatasetFactory::construct(path));
    set_default_names();
}

Dataset::~Dataset() = default;

void Dataset::assign_matrix(Matrix<double>&& m) {
    if (m.M != M) {
        throw except::invalid_operation("Dataset::operator=: Matrix has wrong number of columns. "
        "Expected " + std::to_string(M) + ", but got " + std::to_string(m.M));
    }
    force_assign_matrix(std::move(m));
}

void Dataset::force_assign_matrix(Matrix<double>&& m) {
    this->data = std::move(m.data);
    this->N = m.N;
    this->M = m.M;
}

bool Dataset::empty() const noexcept {
    return data.empty();
}

void Dataset::limit_x(const Limit& limits) {
    if (size() == 0) {return;}
    if (limits.min < x(0) && x(size()-1) < limits.max) {return;}

    Matrix<double> limited(0, M); 
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

    Matrix<double> limited(0, M);
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
    return Matrix::col(index);
}

const ConstColumn<double> Dataset::col(unsigned int index) const {
    return Matrix::col(index);
}

MutableRow<double> Dataset::row(unsigned int index) {
    return Matrix::row(index);
}

const ConstRow<double> Dataset::row(unsigned int index) const {
    return Matrix::row(index);
}

void Dataset::set_col_names(const std::vector<std::string>& names) {
    if (names.size() != M) {throw except::invalid_operation("Dataset::set_col_names: Number of names does not match number of columns. (" + std::to_string(names.size()) + " != " + std::to_string(M) + ")");}
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

void Dataset::save(const io::File& path, const std::string& header) const {
    path.directory().create();

    // check if file was succesfully opened
    std::ofstream output(path);
    if (!output.is_open()) {throw std::ios_base::failure("IntensityFitter::save: Could not open file \"" + path + "\"");}

    // write header
    if (!header.empty()) {
        output << header << std::endl;
    }

    // write column titles
    if (names.size() < M) {throw except::unexpected("Dataset::save: Number of column names (" + std::to_string(names.size()) + ") does not match number of columns (" + std::to_string(M) + ").");}
    for (unsigned int j = 0; j < M; j++) {
        output << std::left << std::setw(16) << names[j] << "\t";
    }
    output << std::endl;

    // write data
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < M-1; j++) {
            output << std::left << std::setw(16) << std::setprecision(8) << std::scientific << index(i, j) << "\t";
        }
        output << index(i, M-1) << "\n";
    }
    output.close();
}

void Dataset::load(const io::ExistingFile& path) {
    auto data = factory::DatasetFactory::construct(path, M);
    if (data->M != M) {throw except::invalid_operation("Dataset::load: Number of columns does not match. (" + std::to_string(data->M) + " != " + std::to_string(M) + ")");}    
    *this = std::move(*data);
    set_default_names();
}

void Dataset::set_default_names() {
    names.resize(M);
    for (unsigned int i = 0; i < M; i++) {
        names[i] = "col_" + std::to_string(i);
    }
}

Dataset Dataset::rolling_average(unsigned int window_size) const {
    Dataset result(*this);
    result.y() = MovingAverage::average_half(y(), window_size);
    return result;
}

Dataset Dataset::interpolate(unsigned int n) const {
    Matrix interpolated(size()*(n+1)-n-1, M);

    std::vector<math::CubicSpline> splines;
    for (unsigned int col_index = 1; col_index < M; ++col_index) {
        splines.push_back(math::CubicSpline(x(), col(col_index)));
    }

    for (unsigned int i = 0; i < size()-1; i++) {
        double x = this->x(i);
        interpolated[i*(n+1)] = row(i);

        double x_next = this->x(i+1);
        double step = (x_next - x)/(n+1);
        for (unsigned int j = 0; j < n; j++) {
            std::vector<double> row_new(M);
            row_new[0] = x + (j+1)*step;;
            for (unsigned int k = 1; k < M; k++) {
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
        return std::vector<double>(M, 0);
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
    Matrix interpolated(newx.size(), M);

    std::vector<math::CubicSpline> splines;
    for (unsigned int col_index = 1; col_index < M; ++col_index) {
        splines.push_back(math::CubicSpline(x(), col(col_index)));
    }

    for (unsigned int i = 0; i < newx.size(); i++) {
        std::vector<double> row_new(M);
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
    if (M != other.M) {throw except::invalid_argument("Dataset::append: Number of columns does not match.");}
    unsigned int n = size();
    extend(other.size());
    for (unsigned int i = 0; i < other.size(); i++) {
        row(n+i) = other.row(i);
    }
}

void Dataset::sort_x() {
    Matrix<double> newdata(N, M);
    std::vector<unsigned int> indices(N);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [this] (unsigned int i, unsigned int j) {return x(i) < x(j);});
    for (unsigned int i = 0; i < N; i++) {
        newdata.row(i) = this->row(indices[i]);
    }
    this->assign_matrix(std::move(newdata));
}

std::string Dataset::to_string() const {
    std::stringstream ss;
    for (unsigned int i = 0; i < size(); i++) {
        for (unsigned int j = 0; j < M; j++) {
            ss << std::setw(16) << std::setprecision(8) << std::scientific << index(i, j) << " ";
        }
        ss << "\n";
    }
    return ss.str();
}

std::vector<unsigned int> Dataset::find_minima(unsigned int min_spacing, double min_prominence) const {
    return math::find_minima(x(), y(), min_spacing, min_prominence);
}

std::vector<unsigned int> Dataset::find_maxima(unsigned int min_spacing, double min_prominence) const {
    return math::find_minima(x(), -y(), min_spacing, min_prominence);
}

unsigned int Dataset::size() const noexcept {
    return size_rows();
}

unsigned int Dataset::size_rows() const noexcept {
    return N;
}

unsigned int Dataset::size_cols() const noexcept {
    return M;
}

bool Dataset::operator==(const Dataset& other) const = default;
Dataset& Dataset::operator=(const Dataset& other) = default;
Dataset& Dataset::operator=(Dataset&& other) = default;