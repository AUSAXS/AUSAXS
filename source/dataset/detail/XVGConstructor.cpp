#include <dataset/detail/XVGConstructor.h>
#include <dataset/Dataset.h>
#include <dataset/SimpleDataset.h>
#include <dataset/Dataset2D.h>
#include <utility/Console.h>
#include <utility/StringUtils.h>
#include <utility/Exceptions.h>
#include <math/Statistics.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>

#include <vector>
#include <string>
#include <fstream>

std::shared_ptr<Dataset> detail::XVGConstructor::construct(const io::ExistingFile& path, unsigned int expected_cols) {
    if (settings::general::verbose) {
        console::print_info("Loading dataset from \"" + path + "\"");
    }

    // check if file was succesfully opened
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("Dataset::load: Could not open file \"" + path + "\"");}

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
        if (skip) {continue;}

        // add values to dataset
        std::vector<double> vals(tokens.size());
        for (unsigned int i = 0; i < tokens.size(); i++) {
            vals[i] = std::stod(tokens[i]);
        }
        row_data.push_back(vals);
        col_number.push_back(vals.size());
    }

    unsigned int mode = stats::mode(col_number);

    // sanity check: comparing the detected number of columns with the header
    if (header.back().find("type") != std::string::npos) {
        std::string type = header.back().substr(6);
        if (mode == 2 && type != "xy") {
            console::print_warning("\tThe column format of the file \"" + path + "\" may be incompatible. Ensure it is of the form [x | y].");
        }
        if (mode == 3 && type != "xydy") {
            console::print_warning("\tThe column format of the file \"" + path + "\" may be incompatible. Ensure it is of the form [x | y | dy].");
        }
        if (mode == 4 && type != "xydxdy") {
            console::print_warning("\tThe column format of the file \"" + path + "\" may be incompatible. Ensure it is of the form [x | y | dy | dx].");
        }
    }

    // check that we have at least the expected number of columns
    if (expected_cols != 0 && mode < expected_cols) {
        throw except::io_error("DATConstructor::construct: File has too few columns. Expected" + std::to_string(expected_cols) + " but found " + std::to_string(mode) + ".");
    }

    // copy the data to the dataset
    std::shared_ptr<Dataset> dataset;
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
        if (mode != expected_cols) {
            // shorten the data to the expected number of columns
            for (unsigned int i = 0; i < data_cols.size(); i++) {
                std::vector<double> row(expected_cols);
                for (unsigned int j = 0; j < expected_cols; j++) {
                    row[j] = data_cols[i][j];
                }
                data_cols[i] = std::move(row);
            }
        }

        // add the data to the dataset
        // dataset = std::make_shared<Dataset>(std::move(data_cols));
        dataset = std::make_shared<Dataset>(0, expected_cols);
        for (unsigned int i = 0; i < data_cols.size(); i++) {
            dataset->push_back(data_cols[i]);
        }
    }

    // skip the first few rows if requested
    if (settings::axes::skip != 0 && settings::general::verbose) {
        std::cout << "\tSkipped " << count - dataset->size() << " data points from beginning of file." << std::endl;
    }

    // verify that at least one row was read correctly
    if (dataset->empty()) {
        throw except::io_error("DATConstructor::construct: No data could be read from the file.");
    }

    // unit conversion
    if (settings::general::verbose) {std::cout << "\tAssuming q is given in units of [nm]." << std::endl;}
    for (unsigned int i = 0; i < dataset->size(); i++) {
        dataset->index(i, 0) /= 10;
    }
    
    // remove all rows outside the specified q-range
    unsigned int N = dataset->size();
    dataset->limit_x(settings::axes::qmin, settings::axes::qmax);
    if (N != dataset->size() && settings::general::verbose) {
        std::cout << "\tRemoved " << N - dataset->size() << " data points outside specified q-range [" << settings::axes::qmin << ", " << settings::axes::qmax << "]." << std::endl;
    }

    // check if the file is abnormally large
    if (dataset->size() > 300) {
        // reread first line
        input.clear();
        input.seekg(0, input.beg);
        getline(input, line);

        // check if file has already been rebinned
        if (line.find("REBINNED") == std::string::npos) {
            // if not, suggest it to the user
            if (settings::general::verbose) {
                std::cout << "\tFile contains more than 300 rows. Consider rebinning the data." << std::endl;
            }
        }
    }

    if (settings::general::verbose) {
        std::cout << "\tSuccessfully read " << dataset->size() << " data points from " << path << std::endl;
    }
    return dataset;
}