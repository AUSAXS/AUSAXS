#include <math/SimpleLeastSquares.h>
#include <utility/Dataset.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <utility/Settings.h>

#include <vector>
#include <string>
#include <fstream>

using std::vector, std::string;

Dataset::Dataset(unsigned int rows, unsigned int cols) : Matrix(rows, cols) {
    for (unsigned int i = 0; i < cols; i++) {
        names.push_back("col" + std::to_string(i));
    }
}

Dataset::Dataset(std::vector<std::vector<double>> cols) : Matrix(cols) {
    for (unsigned int i = 0; i < cols.size(); i++) {
        names.push_back("col " + std::to_string(i));
    }
}

void Dataset::operator=(const Matrix<double>&& other) {
    if (other.M != M) {throw except::invalid_operation("Error in Dataset::operator=: Matrix has wrong number of columns.");}
    this->data = std::move(other.data);
    this->N = other.N;
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
    throw except::invalid_operation("Error in Dataset::col: Column \"" + column + "\" not found. Available columns:\n\t" + utility::join(names, "\n\t"));
}

const ConstColumn<double> Dataset::col(std::string column) const {
    for (size_t i = 0; i < names.size(); ++i) {
        if (names[i] == column) {
            return col(i);
        }
    }
    throw except::invalid_operation("Error in Dataset::col: Column \"" + column + "\" not found. Available columns: " + utility::join(names, "\n"));
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
    if (names.size() != M) {throw except::invalid_operation("Error in Dataset::set_col_names: Number of names does not match number of columns. (" + std::to_string(names.size()) + " != " + std::to_string(M) + ")");}
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
    if (!output.is_open()) {throw std::ios_base::failure("Error in IntensityFitter::save: Could not open file \"" + path + "\"");}

    // write header
    output << header << std::endl;
    for (unsigned int j = 0; j < M; j++) {
        output << std::left << std::setw(14) << names[j] << "\t";
    }
    output << std::endl;

    // write data
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < M-1; j++) {
            output << std::left << std::scientific << std::setw(14) << index(i, j) << "\t";
        }
        output << index(i, M-1) << "\n";
    }
    output.close();
}

#include <math/Statistics.h>
void Dataset::load(std::string path) {
    // check if file was succesfully opened
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("Error in Dataset::load: Could not open file \"" + path + "\"");}

    std::string line;
    std::vector<std::vector<double>> row_data;
    std::vector<unsigned int> col_number;
    while(getline(input, line)) {
        // skip empty lines
        if (line.empty()) {continue;}

        // remove leading whitespace
        vector<string> tokens = utility::split(line, " ,\t\n\r"); // spaces, commas, and tabs can be used as separators

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
            if (tokens[i].find_first_not_of("0123456789+-.Ee\n\r") != string::npos) {
                skip = true;
            }
        }
        if (skip) {continue;}

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
    if (M == 0) {M = mode;}
    else if (M != mode) {throw except::io_error("Error in Dataset::load: Number of columns in file does not match storage capacity of this class. (" + std::to_string(mode) + " != " + std::to_string(M) + ")");}

    // copy all rows with the correct number of columns
    for (unsigned int i = 0; i < row_data.size(); i++) {
        if (row_data[i].size() == mode) {
            push_back(row_data[i]);
        }
    }

    // verify that at least one row was read correctly
    if (size() == 0) {
        throw except::unexpected("Error in Dataset::load: No data was read from the file.");
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
            utility::print_warning("Warning in Dataset::load: File \"" + path + "\" contains more than 300 rows. Consider rebinning it first. (size = " + std::to_string(size()) + ")");
        }
    }    
}