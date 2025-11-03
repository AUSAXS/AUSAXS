// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <dataset/detail/XVGReader.h>
#include <dataset/Dataset.h>
#include <dataset/SimpleDataset.h>
#include <dataset/Dataset2D.h>
#include <utility/Console.h>
#include <utility/StringUtils.h>
#include <utility/Exceptions.h>
#include <math/Statistics.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/Flags.h>

#include <vector>
#include <string>
#include <fstream>

using namespace ausaxs;

std::unique_ptr<Dataset> parse_data(std::vector<std::string>&& header, std::vector<unsigned int>&& col_number, std::vector<std::vector<double>>&& row_data, const io::ExistingFile& path, unsigned int expected_cols) {
    unsigned int mode = stats::mode(col_number);

   // sanity check: comparing the detected number of columns with the header
    if (!header.empty() && header.back().find("type") != std::string::npos) {
        std::string type = header.back().substr(6);
        if (mode == 2) {
            if (type != "xy") {console::print_warning("The column format of the file \"" + path.str() + "\" may be incompatible. Ensure it is of the form [q | I].");}
        }
        else if (mode == 3) {
            if (type != "xydy") {console::print_warning("The column format of the file \"" + path.str() + "\" may be incompatible. Ensure it is of the form [q | I | Ierr].");}
        }
        else if (mode == 4) {
            if (type != "xydxdy") {console::print_warning("The column format of the file \"" + path.str() + "\" may be incompatible. Ensure it is of the form [q | I | Ierr | qerr].");}
        }
    }

    // check that we have at least the expected number of columns
    if (expected_cols != 0 && mode < expected_cols) {
        throw except::io_error("XVGReader::parse_data: File has too few columns. Expected" + std::to_string(expected_cols) + " but found " + std::to_string(mode) + ".");
    }

    // copy the data to the dataset
    std::unique_ptr<Dataset> dataset;
    unsigned int count = 0;
    {
        // first copy all rows with the most common number of columns to a temporary vector
        std::vector<std::vector<double>> data_cols;
        for (unsigned int i = 0; i < row_data.size(); i++) {
            if (row_data[i].size() != mode) {continue;}     // skip rows with the wrong number of columns
            if (count++ < settings::axes::skip) {continue;}  // skip the first few rows if requested
            data_cols.push_back(std::move(row_data[i]));
        }

        // having too many columns is not a problem, but we should inform the user and then ignore the extra columns
        if (expected_cols != 0 && mode != expected_cols) {
            // shorten the data to the expected number of columns
            for (unsigned int i = 0; i < data_cols.size(); i++) {
                std::vector<double> row(expected_cols);
                for (unsigned int j = 0; j < expected_cols; j++) {
                    row[j] = data_cols[i][j];
                }
                data_cols[i] = std::move(row);
            }
            mode = expected_cols;
        }

        // add the data to the dataset
        // dataset = std::make_shared<Dataset>(std::move(data_cols));
        dataset = std::make_unique<Dataset>(0, mode);
        for (unsigned int i = 0; i < data_cols.size(); i++) {
            dataset->push_back(data_cols[i]);
        }
    }

    // skip the first few rows if requested
    if (settings::axes::skip != 0) {
        console::print_text("Skipped " + std::to_string(count - dataset->size_rows()) + " data points from beginning of file.");
    }

    // verify that at least one row was read correctly
    if (dataset->empty()) {
        throw except::io_error("XVGReader::parse_data: No data could be read from the file.");
    }

    // unit conversion
    console::print_text("Assuming q is given in units of [nm].");
    for (unsigned int i = 0; i < dataset->size_rows(); i++) {
        dataset->index(i, 0) /= 10;
    }
    
    // remove all rows outside the specified q-range
    unsigned int N = dataset->size_rows();
    dataset->limit_x(settings::axes::qmin, settings::axes::qmax);
    if (N != dataset->size_rows()) {
        console::print_text(
            "Removed " + std::to_string(N - dataset->size_rows()) + " data points outside specified q-range "
            "[" + std::to_string(settings::axes::qmin) + ", " + std::to_string(settings::axes::qmax) + "]."
        );
    }

    return dataset;
};

std::unique_ptr<Dataset> detail::XVGReader::construct(const io::ExistingFile& path, unsigned int expected_cols) {
    console::print_info("\nReading dataset from \"" + path.str() + "\"");
    console::indent();

    // check if file was succesfully opened
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("XVGReader::parse_data: Could not open file \"" + path.str() + "\"");}

    std::string line;
    std::vector<std::string> header;
    std::vector<std::vector<double>> row_data;
    std::vector<unsigned int> col_number;
    while(getline(input, line)) {
        if (line.empty()) {continue;}                           // skip empty lines
        if (line[0] == '#') {continue;}                         // skip comments
        if (line[0] == '@') {header.push_back(line); continue;} // save header lines

        // remove leading whitespace
        std::vector<std::string> tokens = utility::split(line, " ,\t\n\r"); // spaces, commas, and tabs can be used as separators

        // remove empty tokens
        {
            std::vector<std::string> new_tokens;
            for (unsigned int i = 0; i < tokens.size(); i++) {
                if (tokens[i].empty()) {
                    continue;
                }
                new_tokens.push_back(tokens[i]);
            }
            tokens = std::move(new_tokens);
        }

        // check if all tokens are numbers
        bool skip = false;
        for (unsigned int i = 0; i < tokens.size(); i++) {
            if (tokens[i].find_first_not_of("0123456789+-.Ee\n\r") != std::string::npos) {
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

    auto dataset = parse_data(std::move(header), std::move(col_number), std::move(row_data), path, expected_cols);
    // check if the file is abnormally large
    if (dataset->size_rows() > 300) {
        // reread first line
        input.clear();
        input.seekg(0, input.beg);
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

std::vector<std::unique_ptr<Dataset>> detail::XVGReader::construct_multifile(const io::ExistingFile& path) {
    console::print_info("Reading dataset from \"" + path.str() + "\"");

    // disable text from parsing of each dataset section
    auto verbose = settings::general::verbose;
    settings::general::verbose = false;

    // check if file was succesfully opened
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("XVGReader::construct_multifile: Could not open file \"" + path.str() + "\"");}

    std::string line;
    std::vector<std::string> header;
    std::vector<std::vector<double>> row_data;
    std::vector<unsigned int> col_number;
    std::vector<std::unique_ptr<Dataset>> datasets;
    while(getline(input, line)) {
        if (line.find("@type") == std::string::npos) {continue;}
        else {header.push_back(line);}

        while (getline(input, line) && line[0] != '&') {
            // remove leading whitespace
            std::vector<std::string> tokens = utility::split(line, " ,\t\n\r"); // spaces, commas, and tabs can be used as separators

            // remove empty tokens
            {
                std::vector<std::string> new_tokens;
                for (unsigned int i = 0; i < tokens.size(); i++) {
                    if (tokens[i].empty()) {
                        continue;
                    }
                    new_tokens.push_back(tokens[i]);
                }
                tokens = std::move(new_tokens);
            }

            // check if all tokens are numbers
            bool skip = false;
            for (unsigned int i = 0; i < tokens.size(); i++) {
                if (tokens[i].find_first_not_of("0123456789+-.Ee\n\r") != std::string::npos) {
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

        datasets.push_back(parse_data(std::move(header), std::move(col_number), std::move(row_data), path, 0));
        header.clear();
        row_data.clear();
        col_number.clear();
    }

    settings::general::verbose = verbose;
    console::print_text("Successfully read " + std::to_string(datasets.size()) + " datasets.");
    return datasets;
}