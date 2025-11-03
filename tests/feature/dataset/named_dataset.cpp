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
        std::vector<double> q_values = {0.01, 0.02, 0.03, 0.04, 0.05};
        std::vector<double> intensity = {100.0, 95.0, 85.0, 70.0, 50.0};
        std::vector<double> errors = {2.0, 2.5, 3.0, 3.5, 4.0};
        SimpleDataset raw_data(q_values, intensity, errors);

        // verify internal calculations work with indices
        CHECK(raw_data.x(0) == 0.01);
        CHECK(raw_data.y(0) == 100.0);
        CHECK(raw_data.yerr(0) == 2.0);

        // user-facing output: wrap with NamedDataset for clarity
        NamedSimpleDataset named_output(std::move(raw_data), {"q", "I(q)", "error"});

        // save with meaningful column names (use build directory for portability)
        io::File output_path("temp/tests/named_dataset.dat");
        named_output.save(output_path, "# SAXS experimental data");
        REQUIRE(output_path.exists());

        // read file and verify it has column names
        std::ifstream input(output_path.str());
        std::string header_line;
        std::getline(input, header_line); // comment line
        std::string column_line;
        std::getline(input, column_line); // column names
        input.close();

        CHECK(column_line.find("q") != std::string::npos);
        CHECK(column_line.find("I(q)") != std::string::npos);
        CHECK(column_line.find("error") != std::string::npos);
    }
}

TEST_CASE("NamedDataset: Analysis workflow with named columns") {
    SECTION("Process data with named access") {
        std::vector<double> time = {0.0, 1.0, 2.0, 3.0, 4.0};
        std::vector<double> signal = {10.0, 15.0, 18.0, 16.0, 12.0};
        std::vector<double> baseline = {9.0, 9.1, 9.0, 8.9, 9.0};
        Dataset measurements({time, signal, baseline});
        NamedDataset named_data(std::move(measurements), {"time", "signal", "baseline"});

        // access columns by name for clarity
        auto signal_col = named_data.col("signal");
        auto baseline_col = named_data.col("baseline");

        // calculate corrected signal
        auto corrected = signal_col - baseline_col;
        CHECK_THAT(corrected[0], Catch::Matchers::WithinAbs(1.0, 1e-6));
        CHECK_THAT(corrected[1], Catch::Matchers::WithinAbs(5.9, 1e-6));
        CHECK_THAT(corrected[2], Catch::Matchers::WithinAbs(9.0, 1e-6));
    }
}

TEST_CASE("NamedDataset: Selection by column name") {
    SECTION("Extract specific columns for plotting") {
        // full dataset with many columns
        std::vector<double> q = {0.01, 0.02, 0.03};
        std::vector<double> I_exp = {100.0, 90.0, 75.0};
        std::vector<double> I_err = {2.0, 2.5, 3.0};
        std::vector<double> I_fit = {98.0, 91.0, 76.0};
        std::vector<double> residuals = {2.0, -1.0, -1.0};

        Dataset full_results({q, I_exp, I_err, I_fit, residuals});
        NamedDataset named_results(
            std::move(full_results), 
            {"q", "I_exp", "I_err", "I_fit", "residuals"}
        );
        // select only columns needed for a specific plot
        auto plot_data = named_results.select_columns({"q", "I_exp", "I_fit"});

        CHECK(plot_data.size_cols() == 3);
        CHECK(plot_data.get_col_names() == std::vector<std::string>{"q", "I_exp", "I_fit"});

        // verify data was correctly selected
        CHECK(plot_data.col(0)[0] == 0.01);
        CHECK(plot_data.col(1)[0] == 100.0);
        CHECK(plot_data.col(2)[0] == 98.0);
    }
}

TEST_CASE("NamedDataset: Converting between named and unnamed") {
    SECTION("Start with names, extract data, add different names") {
        Dataset original(3, 2);
        original.col(0) = std::vector<double>{1.0, 2.0, 3.0};
        original.col(1) = std::vector<double>{4.0, 5.0, 6.0};

        NamedDataset input_data(std::move(original), {"input_x", "input_y"});

        std::vector<double> result_col(3);
        for (size_t i = 0; i < 3; ++i) {
            result_col[i] = input_data.col(0)[i]*input_data.col(1)[i];
        }

        // create output with different names
        Dataset output_dataset(3, 2);
        output_dataset.col(0) = input_data.col(0);
        output_dataset.col(1) = result_col;

        NamedDataset output_data(std::move(output_dataset), {"x", "x*y"});

        CHECK(output_data.get_col_names() == std::vector<std::string>{"x", "x*y"});
        CHECK(output_data.col(1)[0] == 4.0);
        CHECK(output_data.col(1)[1] == 10.0);
        CHECK(output_data.col(1)[2] == 18.0);
    }
}