// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/detail/RegularLandscape.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <io/ExistingFile.h>

#include <fstream>

using namespace ausaxs::mini;

void RegularLandscape::rotate90() noexcept {
    x.swap(y);
}

Evaluation RegularLandscape::find_min_val() const {
    if (x.empty() || y.empty()) {
        throw except::size_error("Landscape::find_min_val: x or y is empty.");
    }

    Evaluation min({x[0], y[0]}, z(0, 0));
    for (unsigned int i = 0; i < x.size(); i++) {
        for (unsigned int j = 0; j < y.size(); j++) {
            if (z(i, j) < min.fval) {
                min.fval = z(i, j);
                min.vals = {x[i], y[j]};
            }
        }
    }
    return min;
}

Evaluation RegularLandscape::find_min_eval() const {
    if (evals.empty()) {
        throw except::size_error("Landscape::find_min: Landscape is empty.");
    }

    double min = evals.front().fval;
    auto min_e = evals.front();
    for (const auto& e : evals) {
        if (e.fval < min) {
            min = e.fval;
            min_e = e;
        }
    }
    return min_e;
}

void RegularLandscape::save(std::string filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw except::io_error("Landscape::save: Could not open file " + filename + " for writing.");
    }

    file << "matrix_size: " << x.size() << "\t" << y.size() << "\n";
    for (unsigned int i = 0; i < x.size(); i++) {
        for (unsigned int j = 0; j < y.size(); j++) {
            file << x[i] << "\t" << y[j] << "\t" << z(i, j) << std::endl;
        }
    }

    file.close();
}

RegularLandscape::RegularLandscape(const io::ExistingFile& file) {
    load(file);
}

void RegularLandscape::load(std::string filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw except::io_error("Landscape::load: Could not open file " + filename + " for reading.");
    }

    std::string line;
    std::getline(file, line);
    auto tokens = utility::split(line, " \t");
    if (tokens.size() != 3) {
        throw except::io_error("Landscape::load: Invalid file format.");
    }

    unsigned int x_size = std::stoi(tokens[1]);
    unsigned int y_size = std::stoi(tokens[2]);
    x.resize(x_size);
    y.resize(y_size);
    z = Matrix<double>(x_size, y_size);

    for (unsigned int i = 0; i < x_size; i++) {
        for (unsigned int j = 0; j < y_size; j++) {
            std::getline(file, line);
            tokens = utility::split(line, " \t");
            if (tokens.size() != 3) {
                throw except::io_error("Landscape::load: Invalid file format.");
            }

            // only set each y once
            if (i == 0) {
                y[j] = std::stod(tokens[1]);
            }

            // only set each x once
            if (j == 0) {
                x[i] = std::stod(tokens[0]);
            }

            // always set z
            z(i, j) = std::stod(tokens[2]);
        }
    }

    file.close();
}