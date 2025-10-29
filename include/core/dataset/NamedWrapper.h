// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <dataset/DatasetFwd.h>
#include <math/slices/Slice.h>
#include <io/IOFwd.h>

#include <vector>
#include <string>

namespace ausaxs {
    /**
     * @brief A template wrapper that adds NamedWrapper column functionality to any Dataset class.
     */
    template<typename T>
    struct NamedWrapper {
        T dataset;
        std::vector<std::string> names;

        NamedWrapper() = default;

        template<typename... Args> requires std::is_constructible_v<T, Args...>
        NamedWrapper(Args... args, const std::vector<std::string>& names) : dataset(std::forward<Args>(args)...), names(names) {}
        NamedWrapper(T&& d) : dataset(std::forward<T>(d)) {set_default_names();}
        NamedWrapper(T&& d, const std::vector<std::string>& names) : dataset(std::forward<T>(d)), names(names) {}

        /**
         * @brief Set default column names (col_0, col_1, ...).
         */
        void set_default_names();

        /**
         * @brief Set the column names.
         */
        void set_col_names(const std::vector<std::string>& new_names);

        /**
         * @brief Set a single column name.
         */
        void set_col_names(unsigned int i, const std::string& name);

        /**
         * @brief Get the column names.
         */
        std::vector<std::string> get_col_names() const;

        /**
         * @brief Get a single column name.
         */
        std::string get_col_names(unsigned int i) const;

        /**
         * @brief Check if this dataset has non-default named columns.
         */
        bool is_named() const noexcept;

        /**
         * @brief Get a column based on its name.
         */
        [[nodiscard]] MutableColumn<double> col(std::string_view column);

        /**
         * @brief Get a column based on its name.
         */
        [[nodiscard]] const ConstColumn<double> col(std::string_view column) const;

        /**
         * @brief Create a new dataset with the specified columns.
         */
        NamedWrapper<Dataset> select_columns(const std::vector<std::string>& cols) const;

        /**
         * @brief Get the underlying dataset.
         */
        const T& get_dataset() const noexcept;

        /**
         * @brief Get the underlying dataset.
         */
        T& get_dataset() noexcept;

        /**
         * @brief Write this dataset to the specified file with column names.
         */
        void save(const io::File& path, const std::string& header = "") const;
    };
}