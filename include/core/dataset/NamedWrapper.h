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
    struct NamedWrapper : public T {
        std::vector<std::string> names;

        NamedWrapper() = default;

        NamedWrapper(T&& d) : T(std::forward<T>(d)) {set_default_names();}
        NamedWrapper(T&& d, const std::vector<std::string>& names) : T(std::forward<T>(d)), names(names) {}

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
        using T::col; // bring base class overloads into scope

        /**
         * @brief Create a new dataset with the specified columns.
         */
        NamedWrapper<Dataset> select_columns(std::initializer_list<std::string_view> cols) const;
        using T::select_columns; // bring base class overloads into scope

        /**
         * @brief Write this dataset to the specified file with column names.
         */
        void save(const io::File& path, const std::string& header = "") const;

        static_assert(std::is_base_of_v<Dataset, T>, "NamedWrapper can only be used with Dataset types.");
    };
}