// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dataset/NamedDataset.h>
#include <dataset/SimpleDataset.h>
#include <dataset/Dataset.h>
#include <io/File.h>

#include <fstream>
#include <string>

using namespace ausaxs;

TEST_CASE("NamedDataset: User-facing file I/O") {
    SECTION("Save and load with column names") {
        // Create some experimental data
        std::vector<double> q_values = {0.01, 0.02, 0.03, 0.04, 0.05};
        std::vector<double> intensity = {100.0, 95.0, 85.0, 70.0, 50.0};
        std::vector<double> errors = {2.0, 2.5, 3.0, 3.5, 4.0};
        
        // Internal calculation: use SimpleDataset directly (no names needed)
        SimpleDataset raw_data(q_values, intensity, errors);
        
        // Verify internal calculations work with indices
        CHECK(raw_data.x(0) == 0.01);
        CHECK(raw_data.y(0) == 100.0);
        CHECK(raw_data.yerr(0) == 2.0);
        
        // User-facing output: wrap with NamedDataset for clarity
        NamedDataset<SimpleDataset> named_output(std::move(raw_data), {"q", "I(q)", "error"});
        
        // Save with meaningful column names (use build directory for portability)
        io::File output_path("saxs_data_test.dat");
        
        // Ensure cleanup happens even if test fails
        struct FileGuard {
            io::File path;
            ~FileGuard() { 
                if (path.exists()) {
                    std::remove(path.str().c_str());
                }
            }
        } guard{output_path};
        
        named_output.save(output_path, "# SAXS experimental data");
        
        // Verify file was created
        REQUIRE(output_path.exists());
        
        // Read file and verify it has column names
        std::ifstream input(output_path.str());
        std::string header_line;
        std::getline(input, header_line); // Comment line
        std::string column_line;
        std::getline(input, column_line); // Column names
        input.close();
        
        CHECK(column_line.find("q") != std::string::npos);
        CHECK(column_line.find("I(q)") != std::string::npos);
        CHECK(column_line.find("error") != std::string::npos);
        
        // Cleanup handled by FileGuard destructor
    }
}

TEST_CASE("NamedDataset: Analysis workflow with named columns") {
    SECTION("Process data with named access") {
        // Start with a dataset that has meaningful column names
        std::vector<double> time = {0.0, 1.0, 2.0, 3.0, 4.0};
        std::vector<double> signal = {10.0, 15.0, 18.0, 16.0, 12.0};
        std::vector<double> baseline = {9.0, 9.1, 9.0, 8.9, 9.0};
        
        Dataset measurements({time, signal, baseline});
        NamedDataset<Dataset> named_data(std::move(measurements), {"time", "signal", "baseline"});
        
        // Access columns by name for clarity
        auto signal_col = named_data.col("signal");
        auto baseline_col = named_data.col("baseline");
        
        // Calculate corrected signal
        std::vector<double> corrected(5);
        for (size_t i = 0; i < 5; ++i) {
            corrected[i] = signal_col[i] - baseline_col[i];
        }
        
        CHECK_THAT(corrected[0], Catch::Matchers::WithinRel(1.0, 0.01));
        CHECK_THAT(corrected[1], Catch::Matchers::WithinRel(5.9, 0.01));
        CHECK_THAT(corrected[2], Catch::Matchers::WithinRel(9.0, 0.01));
    }
}

TEST_CASE("NamedDataset: Selection by column name") {
    SECTION("Extract specific columns for plotting") {
        // Full dataset with many columns
        std::vector<double> q = {0.01, 0.02, 0.03};
        std::vector<double> I_exp = {100.0, 90.0, 75.0};
        std::vector<double> I_err = {2.0, 2.5, 3.0};
        std::vector<double> I_fit = {98.0, 91.0, 76.0};
        std::vector<double> residuals = {2.0, -1.0, -1.0};
        
        Dataset full_results({q, I_exp, I_err, I_fit, residuals});
        NamedDataset<Dataset> named_results(std::move(full_results), 
                                            {"q", "I_exp", "I_err", "I_fit", "residuals"});
        
        // Select only columns needed for a specific plot
        auto plot_data = named_results.select_columns({"q", "I_exp", "I_fit"});
        
        CHECK(plot_data.dataset.size_cols() == 3);
        CHECK(plot_data.get_col_names() == std::vector<std::string>{"q", "I_exp", "I_fit"});
        
        // Verify data was correctly selected
        CHECK(plot_data.dataset.col(0)[0] == 0.01);
        CHECK(plot_data.dataset.col(1)[0] == 100.0);
        CHECK(plot_data.dataset.col(2)[0] == 98.0);
    }
}

TEST_CASE("NamedDataset: Converting between named and unnamed") {
    SECTION("Start with names, extract data, add different names") {
        // User provides named data
        Dataset original(3, 2);
        original.col(0) = std::vector<double>{1.0, 2.0, 3.0};
        original.col(1) = std::vector<double>{4.0, 5.0, 6.0};
        
        NamedDataset<Dataset> input_data(std::move(original), {"input_x", "input_y"});
        
        // Perform internal calculation (use dataset directly, names not needed)
        auto& calc_data = input_data.dataset;
        std::vector<double> result_col(3);
        for (size_t i = 0; i < 3; ++i) {
            result_col[i] = calc_data.col(0)[i] * calc_data.col(1)[i];
        }
        
        // Create output with different names
        Dataset output_dataset(3, 2);
        output_dataset.col(0) = calc_data.col(0);
        output_dataset.col(1) = result_col;
        
        NamedDataset<Dataset> output_data(std::move(output_dataset), {"x", "x*y"});
        
        CHECK(output_data.get_col_names() == std::vector<std::string>{"x", "x*y"});
        CHECK(output_data.dataset.col(1)[0] == 4.0);
        CHECK(output_data.dataset.col(1)[1] == 10.0);
        CHECK(output_data.dataset.col(1)[2] == 18.0);
    }
}
