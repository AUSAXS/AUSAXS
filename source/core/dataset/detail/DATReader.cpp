// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <dataset/detail/DATReader.h>

#include <dataset/Dataset.h>
#include <math/Statistics.h>
#include <settings/Flags.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <utility/Console.h>
#include <utility/StringUtils.h>

#include <fstream>
#include <string>
#include <vector>

using namespace ausaxs;

std::unique_ptr<Dataset> detail::DATReader::construct(const io::ExistingFile& path, int expected_cols) {
    console::print_info("\nReading dataset from \"" + path.str() + "\"");
    console::indent();

    // check if file was succesfully opened
    std::ifstream input(path);
    if (!input.is_open()) {throw ausaxs::except::io_error("DATReader::construct: Could not open file \"" + path.str() + "\"");}

    std::string line;
    std::vector<std::string> header;
    std::vector<std::vector<double>> row_data;
    std::vector<int> col_number;
    while(getline(input, line)) {
        // skip empty lines
        if (line.empty()) {continue;}

        // remove leading whitespace
        std::vector<std::string> tokens = utility::split(line, " ,\t\n\r"); // spaces, commas, and tabs can be used as separators

        // remove empty tokens
        {
            std::vector<std::string> new_tokens;
            for (const auto& token : tokens) {
                if (token.empty()) {
                    continue;
                }
                new_tokens.push_back(token);
            }
            tokens = std::move(new_tokens);
        }

        // check if all tokens are numbers
        bool skip = false;
        for (const auto& token : tokens) {
            if (token.find_first_not_of("0123456789+-.Ee\n\r") != std::string::npos) {
                skip = true;
            }
        }
        if (skip) {
            header.push_back(line);
            continue;
        }

        // add values to dataset
        std::vector<double> vals(tokens.size());
        for (int i = 0; i < static_cast<int>(tokens.size()); i++) {
            vals[i] = std::stod(tokens[i]);
        }
        row_data.push_back(vals);
        col_number.push_back(static_cast<int>(vals.size()));
    }

    // determine the most common number of columns, since that will likely be the data
    int mode = stats::mode(col_number);
    switch (mode) {
        case 2: 
            console::print_text("2 columns detected. Assuming the format is [q | I]");
            break;
        case 3:
            console::print_text("3 columns detected. Assuming the format is [q | I | Ierr]");
            break;
        default:
            console::print_text(std::to_string(mode) + " columns detected. Assuming the format is [q | I | Ierr | ... ]");
    }

    // check that we have at least the expected number of columns
    if (expected_cols != 0 && mode < expected_cols) {
        throw except::io_error("DATReader::construct: File has too few columns. Expected" + std::to_string(expected_cols) + " but found " + std::to_string(mode) + ".");
    }

    // copy the data to the dataset
    std::unique_ptr<Dataset> dataset;
    int count = 0;
    {
        // first copy all rows with the most common number of columns to a temporary vector
        std::vector<std::vector<double>> data_cols;
        for (auto& row : row_data) {
            if (static_cast<int>(row.size()) != mode) {continue;}     // skip rows with the wrong number of columns
            if (count++ < settings::axes::skip) {continue;}  // skip the first few rows if requested
            data_cols.emplace_back(std::move(row));
        }

        // having too many columns is not a problem, but we should inform the user and then ignore the extra columns
        if (expected_cols != 0 && mode != expected_cols) {
            // shorten the data to the expected number of columns
            for (auto& col : data_cols) {
                std::vector<double> row(expected_cols);
                for (int j = 0; j < expected_cols; j++) {
                    row[j] = col[j];
                }
                col = std::move(row);
            }
            mode = expected_cols; // update mode to the number of columns we actually have
        }

        // add the data to the dataset
        // dataset = std::make_shared<Dataset>(std::move(data_cols));
        dataset = std::make_unique<Dataset>(0, mode);
        for (const auto& col : data_cols) {
            dataset->push_back(col);
        }
    }

    // skip the first few rows if requested
    if (settings::axes::skip != 0) {
        console::print_text("Skipped " + std::to_string(count - dataset->size_rows()) + " data points from beginning of file.");
    }

    // verify that at least one row was read correctly
    if (dataset->empty()) {
        throw except::io_error("DATReader::construct: No data could be read from the file.");
    }

    // scan the headers for units. must be either [Å] or [nm]
    settings::general::QUnit unit = settings::general::input_q_unit;
    if (!settings::general::helper::is_user_defined(unit)) {
        bool found = false;
        for (auto& s : header) {
            if (s.starts_with("[nm]") || s.starts_with("[nm^-1]")) {
                if (!settings::general::helper::is_nanometers(settings::general::input_q_unit)) {
                    console::print_warning("Warning: File contains unit [nm], but default is set to [A]. Assuming [nm] is correct.");
                }
                unit = settings::general::QUnit::NM;
                found = true;
                break;
            }
            if ((s.starts_with("[A]") || s.starts_with("[A^-1]"))) {
                if (settings::general::helper::is_nanometers(settings::general::input_q_unit)) {
                    console::print_warning("Warning: File contains unit [A], but default is set to [nm]. Assuming [A] is correct.");
                }
                unit = settings::general::QUnit::A;
                found = true;
                break;
            }
        }
        if (!found) {
            if (1 < dataset->x().back()) {
                console::print_text("Detected q-values larger than 1. Assuming the unit is [nm].");
                unit = settings::general::QUnit::NM;
            } else if (dataset->x().front() < 1e-3) {
                console::print_text("Detected q-values smaller than 1e-3. Assuming the unit is [A].");
                unit = settings::general::QUnit::A;
            } else {
                console::print_warning("Warning: The q-value unit is ambiguous. Assuming [A].");
            }
        }
    }
    if (settings::general::helper::is_nanometers(unit)) {
        console::print_text("Scaling all q-values by 1/10 to convert from inverse [nm] to inverse [A].");
        for (int i = 0; i < dataset->size_rows(); i++) {
            dataset->index(i, 0) /= 10;
        }
        settings::flags::last_parsed_unit = static_cast<char>(settings::general::QUnit::NM);
    }

    // remove all rows outside the specified q-range
    if (settings::axes::clamp_to_qrange) {
        int N = dataset->size_rows();
        dataset->limit_x(settings::axes::qmin, settings::axes::qmax);
        if (N != dataset->size_rows()) {
            console::print_text(
                "Removed " + std::to_string(N - dataset->size_rows()) + " data points outside specified q-range "
                "[" + std::to_string(settings::axes::qmin) + ", " + std::to_string(settings::axes::qmax) + "]."
            );
        }
    }

    // check if the file is abnormally large
    if (dataset->size_rows() > 300) {
        // reread first line
        input.clear();
        input.seekg(0, std::ifstream::beg);
        getline(input, line);

        // check if file has already been rebinned
        if (!settings::flags::data_rebin) {
            // if not, suggest it to the user
            console::print_text_minor("File contains more than 300 rows. Consider rebinning the data.");
        }
    }

    console::print_text("Successfully read " + std::to_string(dataset->size_rows()) + " data points");
    console::unindent();
    return dataset;
}