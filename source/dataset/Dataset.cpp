#include <math/SimpleLeastSquares.h>
#include <math/Statistics.h>
#include <math/MovingAverager.h>
#include <math/CubicSpline.h>
#include <dataset/Dataset.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <dataset/DatasetFactory.h>

#include <vector>
#include <string>
#include <fstream>

Dataset::Dataset(unsigned int rows, unsigned int cols) : Matrix(rows, cols) {
    set_default_names();
}

Dataset::Dataset(std::vector<std::vector<double>> cols) : Matrix(cols) {
    set_default_names();
}

Dataset::Dataset(std::string path) : Dataset() {
    *this = std::move(*factory::DatasetFactory::construct(path));
    set_default_names();
}

void Dataset::assign_matrix(const Matrix<double>&& m) {
    if (m.M != M) {
        throw except::invalid_operation("Dataset::operator=: Matrix has wrong number of columns. "
        "Expected " + std::to_string(M) + ", but got " + std::to_string(m.M));
    }
    force_assign_matrix(std::move(m));
}

void Dataset::force_assign_matrix(const Matrix<double>&& m) {
    this->data = std::move(m.data);
    this->N = m.N;
    this->M = m.M;
}

std::size_t Dataset::size() const noexcept {
    return N;
}

bool Dataset::empty() const noexcept {
    return size() == 0;
}

void Dataset::limit_x(const Limit& limits) {
    if (size() == 0) {return;}
    if (limits.min < x(0) && x(size()-1) < limits.max) {return;}

    Matrix<double> limited(0, M); 
    for (unsigned int i = 0; i < size(); i++) {
        double val = x(i);
        if (val < limits.min || limits.max < val) {continue;}
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

Column<double> Dataset::col(std::string column) {
    for (size_t i = 0; i < names.size(); ++i) {
        if (names[i] == column) {
            return col(i);
        }
    }
    throw except::invalid_operation("Dataset::col: Column \"" + column + "\" not found. Available columns:\n\t" + utility::join(names, "\n\t"));
}

const ConstColumn<double> Dataset::col(std::string column) const {
    for (size_t i = 0; i < names.size(); ++i) {
        if (names[i] == column) {
            return col(i);
        }
    }
    throw except::invalid_operation("Dataset::col: Column \"" + column + "\" not found. Available columns: " + utility::join(names, "\n"));
}

Column<double> Dataset::col(unsigned int index) {
    return Matrix::col(index);
}

const ConstColumn<double> Dataset::col(unsigned int index) const {
    return Matrix::col(index);
}

Row<double> Dataset::row(unsigned int index) {
    return Matrix::row(index);
}

const ConstRow<double> Dataset::row(unsigned int index) const {
    return Matrix::row(index);
}

void Dataset::set_col_names(std::vector<std::string> names) {
    if (names.size() != M) {throw except::invalid_operation("Dataset::set_col_names: Number of names does not match number of columns. (" + std::to_string(names.size()) + " != " + std::to_string(M) + ")");}
    this->names = names;
}

void Dataset::set_col_names(unsigned int i, std::string name) {
    names[i] = name;
}

std::vector<std::string> Dataset::get_col_names() {
    return names;
}

std::string Dataset::get_col_names(unsigned int i) {
    return names[i];
}

void Dataset::save(std::string path, std::string header) const {
    utility::create_directory(path);

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

void Dataset::load(std::string path) {
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

void Dataset::interpolate(unsigned int n) {
    CubicSpline spline(x().to_vector(), y().to_vector());
    Matrix interpolated(size()*(n+1)-n-1, 2);
    for (unsigned int i = 0; i < size()-1; i++) {
        double x = this->x(i);
        double y = this->y(i);
        interpolated[i*(n+1)] = {x, y}; 

        double x_next = this->x(i+1);
        double step = (x_next - x)/(n+1);
        for (unsigned int j = 0; j < n; j++) {
            double x_new = x + (j+1)*step;
            double y_new = spline.spline(x_new);
            interpolated[i*(n+1) + j + 1] = {x_new, y_new};
        }
    }
    // force assign since we can only interpolate the x and y columns. 
    force_assign_matrix(std::move(interpolated));
}

void Dataset::append(const Dataset& other) {
    if (M != other.M) {throw except::invalid_argument("Dataset::append: Number of columns does not match.");}
    unsigned int n = size();
    extend(other.size());
    for (unsigned int i = 0; i < other.size(); i++) {
        row(n+i) = other.row(i);
    }
}

#include <sstream>
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