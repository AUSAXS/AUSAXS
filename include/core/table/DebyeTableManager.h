// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <table/DebyeTable.h>
#include <constants/ConstantsAxes.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <type_traits>

namespace ausaxs::table {
    class DebyeTableManager {
        public:
            DebyeTableManager();
            DebyeTableManager(const DebyeTableManager&);
            DebyeTableManager(DebyeTableManager&&) noexcept;
            DebyeTableManager& operator=(const DebyeTableManager&);
            DebyeTableManager& operator=(DebyeTableManager&&) noexcept;

            /**
             * @brief Get the sinc(x) lookup table for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table() const;

            void reset_to_default();

            template<typename T, typename = std::enable_if_t<std::disjunction_v<
                std::is_rvalue_reference<T&&>,
                std::is_same<T, const std::vector<double>&>,
                std::is_same<T, std::vector<double>&>
            >>>
            void set_q_axis(T&& q_axis);

            template<typename T, typename = std::enable_if_t<std::disjunction_v<
                std::is_rvalue_reference<T&&>,
                std::is_same<T, const std::vector<double>&>,
                std::is_same<T, std::vector<double>&>
            >>>
            void set_d_axis(T&& d_axis);

        private:
            struct {
                std::vector<double> axis;
                bool defaulted = true;
            } q;
            struct {
                std::vector<double> axis;
                bool defaulted = true;
            } d;
            mutable std::unique_ptr<table::DebyeTable> custom_sinc_table;
            mutable bool recalculate = true; // whether the custom table needs to be recalculated
            bool use_custom_table = false;
    };
}