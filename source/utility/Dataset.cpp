#include <math/SimpleLeastSquares.h>
#include <math/Statistics.h>
#include <math/MovingAverager.h>
#include <math/CubicSpline.h>
#include <utility/Dataset.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <utility/Settings.h>

#include <vector>
#include <string>
#include <fstream>

Dataset::Dataset(unsigned int rows, unsigned int cols) : Matrix(rows, cols) {
    set_default_names();
}

Dataset::Dataset(std::vector<std::vector<double>> cols) : Matrix(cols) {
    set_default_names();
}

void Dataset::assign_matrix(const Matrix<double>&& m) {
    if (m.M != M) {throw except::invalid_operation("Dataset::operator=: Matrix has wrong number of columns.");}
    this->data = std::move(m.data);
    this->N = m.N;
}

std::size_t Dataset::size() const noexcept {
    return N;
}

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
    if (setting::general::verbose) {
        utility::print_info("Loading dataset from \"" + path + "\"");
    }

    // check if file was succesfully opened
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("Dataset::load: Could not open file \"" + path + "\"");}

    std::string line;
    std::vector<std::string> header;
    std::vector<std::vector<double>> row_data;
    std::vector<unsigned int> col_number;
    while(getline(input, line)) {
        // skip empty lines
        if (line.empty()) {continue;}

        // remove leading whitespace
        std::vector<std::string> tokens = utility::split(line, " ,\t\n\r"); // spaces, commas, and tabs can be used as separators

        // remove empty tokens
        for (unsigned int i = 0; i < tokens.size(); i++) {
            if (tokens[i].empty()) {
                tokens.erase(tokens.begin() + i);
                i--;
            }
        }

        // check if all tokens are numbers
        bool skip = false;
        for (unsigned int i = 0; i < tokens.size(); i++) {
            if (tokens[i].find_first_not_of("0123456789+-.Ee\n\r") != std::string::npos) {
                skip = true;
            }
        }
        if (skip) {
            header.push_back(line);
            continue;
        }

        // add values to dataset
        std::vector<double> vals(tokens.size());
        for (unsigned int i = 0; i < tokens.size(); i++) {
            vals[i] = std::stod(tokens[i]);
        }
        row_data.push_back(vals);
        col_number.push_back(vals.size());
    }

    // determine the most common number of columns, since that will likely be the data
    unsigned int mode = stats::mode(col_number);
    if (M == 0) {
        M = mode;
        set_default_names();
    }
    else if (M < mode) {throw except::io_error("Dataset::load: Number of columns in file does not match storage capacity of this class. (" + std::to_string(mode) + " != " + std::to_string(M) + ")");}

    // check if the file has the same number of columns as the dataset
    if (mode < M) {
        // if the file has less columns, fill the remaining columns with zeros
        for (unsigned int i = 0; i < M-mode; i++) {
            for (unsigned int i = 0; i < row_data.size(); i++) {
                row_data[i].push_back(0);
            }
        }
        mode = M;
    }

    // copy all rows with the correct number of columns
    unsigned int count = 0;
    for (unsigned int i = 0; i < row_data.size(); i++) {
        if (row_data[i].size() != mode) {continue;}
        if (count++ < setting::axes::skip) {continue;}
        push_back(row_data[i]);
    }
    if (setting::axes::skip != 0 && setting::general::verbose) {
        std::cout << "\tSkipped " << count - size() << " data points from beginning of file." << std::endl;
    }

    // verify that at least one row was read correctly
    if (size() == 0) {
        throw except::unexpected("Dataset::load: No data could be read from the file.");
    }

    // scan the headers for units. must be either [Å] or [nm]
    bool found_unit = false;
    for (auto& s : header) {
        if (s.find("[nm]") != std::string::npos) {
            std::cout << "\tUnit [nm] detected. Scaling all q values by 1/10." << std::endl;
            for (unsigned int i = 0; i < size(); i++) {
                index(i, 0) /= 10;
            }
            found_unit = true;
        } else if (s.find("[Å]") != std::string::npos) {
            std::cout << "\tUnit [Å] detected. No scaling necessary." << std::endl;
            found_unit = true;
        }
    }
    if (!found_unit) {
        std::cout << "\tNo unit detected. Assuming [Å]." << std::endl;
    }

    // check if the file is abnormally large
    if (size() > 300) {
        // reread first line
        input.clear();
        input.seekg(0, input.beg);
        getline(input, line);

        // check if file has already been rebinned
        if (line.find("REBINNED") == std::string::npos) {
            // if not, suggest it to the user
            if (setting::general::verbose) {
                std::cout << "\tFile contains more than 300 rows. Consider rebinning the data." << std::endl;
            }
        }
    }

    if (setting::general::verbose) {
        std::cout << "\tSuccessfully read " << size() << " data points from " << path << std::endl;
    }
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
    *this = Dataset(std::move(interpolated));
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