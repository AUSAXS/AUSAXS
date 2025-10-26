// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Matrix.h>
#include <io/File.h>
#include <utility/Exceptions.h>

#include <vector>
#include <string>
#include <utility>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace ausaxs {
    /**
     * @brief A template wrapper that adds named column functionality to any Dataset class.
     * 
     * This wrapper allows any dataset class to have named columns without requiring
     * the base dataset to always store column names. This is useful for user-facing
     * operations where column names provide clarity (e.g., saving files, plotting),
     * while internal operations can use the unwrapped dataset for efficiency.
     * 
     * Usage patterns:
     * - Use NamedDataset<T> for user-facing code where column names add clarity
     *   (e.g., file I/O, plotting, reporting results)
     * - Use T directly (Dataset, SimpleDataset, etc.) for internal calculations
     *   where column indices are sufficient
     * 
     * Example:
     * @code
     *   // Internal calculation - use Dataset directly
     *   Dataset raw_data(rows, cols);
     *   raw_data.y() = some_calculation(raw_data.x());
     * 
     *   // User-facing output - wrap with NamedDataset
     *   NamedDataset<Dataset> named_data(std::move(raw_data), {"q", "I", "I_err"});
     *   named_data.save("output.dat");
     * @endcode
     * 
     * @tparam T The dataset type to wrap (e.g., Dataset, SimpleDataset, Dataset2D).
     */
    template<typename T>
    struct NamedDataset {
        T dataset;
        std::vector<std::string> names;

        /**
         * @brief Default constructor.
         */
        NamedDataset() = default;

        /**
         * @brief Construct a NamedDataset from a dataset with default column names.
         * 
         * @param d The dataset to wrap.
         */
        NamedDataset(T&& d) : dataset(std::forward<T>(d)) {
            set_default_names();
        }

        /**
         * @brief Construct a NamedDataset from a dataset with explicit column names.
         * 
         * @param d The dataset to wrap.
         * @param names The column names.
         */
        NamedDataset(T&& d, const std::vector<std::string>& names) : dataset(std::forward<T>(d)), names(names) {}

        /**
         * @brief Construct a NamedDataset from a dataset (copy) with default column names.
         * 
         * @param d The dataset to wrap.
         */
        NamedDataset(const T& d) : dataset(d) {
            set_default_names();
        }

        /**
         * @brief Construct a NamedDataset from a dataset (copy) with explicit column names.
         * 
         * @param d The dataset to wrap.
         * @param names The column names.
         */
        NamedDataset(const T& d, const std::vector<std::string>& names) : dataset(d), names(names) {}

        /**
         * @brief Set default column names (col_0, col_1, ...).
         */
        void set_default_names() {
            names.resize(dataset.size_cols());
            for (unsigned int i = 0; i < dataset.size_cols(); i++) {
                names[i] = "col_" + std::to_string(i);
            }
        }

        /**
         * @brief Set the column names.
         * 
         * @param new_names The new column names.
         * @throw except::invalid_operation if the number of names doesn't match the number of columns.
         */
        void set_col_names(const std::vector<std::string>& new_names) {
            if (new_names.size() != dataset.size_cols()) {
                throw except::invalid_operation(
                    "NamedDataset::set_col_names: Number of names does not match number of columns. "
                    "(" + std::to_string(new_names.size()) + " != " + std::to_string(dataset.size_cols()) + ")"
                );
            }
            names = new_names;
        }

        /**
         * @brief Set a single column name.
         * 
         * @param i The column index.
         * @param name The new name.
         */
        void set_col_names(unsigned int i, const std::string& name) {
            names[i] = name;
        }

        /**
         * @brief Get the column names.
         * 
         * @return The column names.
         */
        std::vector<std::string> get_col_names() const {
            return names;
        }

        /**
         * @brief Get a single column name.
         * 
         * @param i The column index.
         * @return The column name.
         */
        std::string get_col_names(unsigned int i) const {
            return names[i];
        }

        /**
         * @brief Check if this dataset has non-default named columns.
         * 
         * @return true if any column has a name other than the default "col_N" format.
         */
        bool is_named() const noexcept {
            for (unsigned int i = 0; i < dataset.size_cols(); i++) {
                if (names[i] != "col_" + std::to_string(i)) {
                    return true;
                }
            }
            return false;
        }

        /**
         * @brief Get a column based on its name.
         */
        [[nodiscard]] MutableColumn<double> col(std::string_view column) {
            for (unsigned int i = 0; i < names.size(); ++i) {
                if (names[i] == column) {
                    return dataset.col(i);
                }
            }
            throw except::invalid_argument("NamedDataset::col: Column \"" + std::string(column) + "\" not found.");
        }

        /**
         * @brief Get a column based on its name.
         */
        [[nodiscard]] const ConstColumn<double> col(std::string_view column) const {
            for (unsigned int i = 0; i < names.size(); ++i) {
                if (names[i] == column) {
                    return dataset.col(i);
                }
            }
            throw except::invalid_argument("NamedDataset::col: Column \"" + std::string(column) + "\" not found.");
        }

        /**
         * @brief Create a new dataset with the specified columns.
         */
        NamedDataset<T> select_columns(const std::vector<std::string>& cols) const {
            std::vector<unsigned int> indices(cols.size());
            std::transform(cols.begin(), cols.end(), indices.begin(), [this](const std::string& name) {
                for (unsigned int i = 0; i < names.size(); i++) {
                    if (names[i] == name) {return i;}
                }
                throw except::invalid_argument("NamedDataset::select_columns: Column \"" + name + "\" not found.");
            });
            
            T new_dataset = dataset.select_columns(indices);
            std::vector<std::string> col_names(cols.size());
            for (unsigned int i = 0; i < cols.size(); i++) {
                col_names[i] = names[indices[i]];
            }
            return NamedDataset<T>(std::move(new_dataset), col_names);
        }

        /**
         * @brief Write this dataset to the specified file with column names.
         * 
         * @param path The path to the save location.
         * @param header The header for the file.
         */
        void save(const io::File& path, const std::string& header = "") const {
            path.directory().create();

            std::ofstream output(path);
            if (!output.is_open()) {
                throw std::ios_base::failure("NamedDataset::save: Could not open file \"" + path.str() + "\"");
            }

            // write header
            if (!header.empty()) {
                output << header << std::endl;
            }

            // write column titles
            if (names.size() < dataset.size_cols()) {
                throw except::unexpected(
                    "NamedDataset::save: Number of column names (" + std::to_string(names.size()) + ") "
                    "does not match number of columns (" + std::to_string(dataset.size_cols()) + ")."
                );
            }
            for (unsigned int j = 0; j < dataset.size_cols(); j++) {
                output << std::left << std::setw(16) << names[j] << "\t";
            }
            output << std::endl;

            // write data
            for (unsigned int i = 0; i < dataset.size_rows(); i++) {
                for (unsigned int j = 0; j < dataset.size_cols()-1; j++) {
                    output << std::left << std::setw(16) << std::setprecision(8) << std::scientific << dataset.index(i, j) << "\t";
                }
                output << dataset.index(i, dataset.size_cols()-1) << "\n";
            }
            output.close();
        }
    };
}
