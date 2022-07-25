#include <math/SimpleLeastSquares.h>
#include <utility/Dataset.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <utility/Settings.h>

#include <vector>
#include <string>
#include <fstream>

using std::vector, std::string;

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
    throw except::invalid_operation("Error in Dataset::col: Column not found.");
}

const ConstColumn<double> Dataset::col(std::string column) const {
    for (size_t i = 0; i < names.size(); ++i) {
        if (names[i] == column) {
            return col(i);
        }
    }
    throw except::invalid_operation("Error in Dataset::col: Column not found.");
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

void Dataset::save(std::string path) const {
    utility::create_directories(path);

    // check if file was succesfully opened
    std::ofstream output(path);
    if (!output.is_open()) {throw std::ios_base::failure("Error in IntensityFitter::save: Could not open file \"" + path + "\"");}

    // prepare header & writer function
    std::function<string(unsigned int)> writer;
    string header = "broken header\n";
    if (M == 2) {
        header = "x y\n";
        writer = [this] (unsigned int i) {return std::to_string(index(i, 0)) + " " + std::to_string(index(i, 1)) + "\n";};
    } else if (M == 3) {
        header = "x y yerr\n";
        writer = [this] (unsigned int i) {return std::to_string(index(i, 0)) + " " + std::to_string(index(i, 1)) + " " + std::to_string(index(i, 2)) + "\n";};
    } else if (M == 4) {
        header = "x y yerr xerr\n";
        writer = [this] (unsigned int i) {return std::to_string(index(i, 0)) + " " + std::to_string(index(i, 1)) + " " + std::to_string(index(i, 2)) + " " + std::to_string(index(i, 3)) + "\n";};
    } else {
        throw except::invalid_operation("Error in IDataset::save: Dataset has wrong number of columns.");
    }

    // write to disk
    output << header;
    for (unsigned int i = 0; i < size(); i++) {
        output << writer(i);
    }
    output.close();
}

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
void Dataset::load(std::string path) {
    // check if file was succesfully opened
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("Error in IDataset::load: Could not open file \"" + path + "\"");}

    std::string line; // placeholder for the current line
    while(getline(input, line)) {
        if (line[0] == ' ') {line = line.substr(1);} // fix leading space
        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of(" ,\t")); // spaces, commas, and tabs can all be used as separators (but not a mix of them)

        // remove empty tokens
        for (unsigned int i = 0; i < tokens.size(); i++) {
            if (tokens[i].empty()) {tokens.erase(tokens.begin()+i);}
        }

        // determine if we are in some sort of header
        if (tokens.size() < 2 || tokens.size() > 4) {continue;} // too many separators
        bool skip = false;
        for (unsigned int i = 0; i < tokens.size(); i++) { // check if all tokens are numbers
            if (tokens[i].find_first_not_of("0123456789-.Ee\n\r") != string::npos) {
                skip = true;
            }
        }
        if (skip) {continue;}

        // now we are most likely beyond any headers
        double _q, _I, _sigma;
        _q = std::stod(tokens[0]); // we know for sure that the strings are convertible to numbers (boost check)
        _I = std::stod(tokens[1]);
        if (_q > 10) {continue;} // probably not a q-value if it's larger than 10

        // check user-defined limits
        if (_q < setting::fit::q_low) {continue;}
        if (_q > setting::fit::q_high) {continue;}

        // add the values to our vectors
        // this is a fair bit more complicated than strictly necessary
        if (tokens.size() == 4) {
            if (M == 4) {
                push_back({_q, _I, std::stod(tokens[2]), std::stod(tokens[3])});
            } else if (M == 3) {
                throw except::invalid_operation("Error in IDataset::load: File has four columns, but a SimpleDataset only supports three. Use a Dataset instance instead.");
            } else {
                throw except::unexpected("Error in IDataset::load: Unknown data layout.");
            }
        }
        else if (tokens.size() == 3) {
            if (M == 4) {
                push_back({_q, _I, std::stod(tokens[2]), 0});
            } else if (M == 3) {
                push_back({_q, _I, std::stod(tokens[2])});
            } else {
                throw except::unexpected("Error in IDataset::load: Unknown data layout.");
            }
        }
        else if (tokens.size() == 2) {
            if (M == 4) {
                push_back({_q, _I, 0, 0});
            } else if (M == 3) {
                push_back({_q, _I, 0});
            } else {
                throw except::unexpected("Error in IDataset::load: Unknown data layout.");
            }
        } else {
            throw except::unexpected("Error in IDataset::load: Unknown data layout.");
        }
    }
    input.close();
    if (size() == 0) {throw except::unexpected("Error in IDataset::load: No data was read from the file.");}
}