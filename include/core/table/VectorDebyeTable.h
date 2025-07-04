// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <table/Table.h>
#include <table/DebyeTable.h>
#include <utility/Concepts.h>

namespace ausaxs::table {
    class VectorDebyeTable : public DebyeTable, private Table {
        public:
            VectorDebyeTable();

            /**
             * @brief Initialize a new runtime sinc lookup table for the given d-axis. 
             *        The default q-axis from costants::axes::q_axis will be used.
             */
            VectorDebyeTable(const std::vector<constants::axes::d_type>& d);

            /**
             * @brief Initialize a new runtime sinc lookup table for the given d-axis. 
             *        The default q-axis from costants::axes::q_axis will be used.
             */
            VectorDebyeTable(const std::array<constants::axes::d_type, constants::axes::d_axis.bins>& d);

            /**
             * @brief Initialize a new runtime sinc lookup table for the given d-axis and q-axis. 
             */
            VectorDebyeTable(const std::vector<constants::axes::d_type>& d, const std::vector<double>& q);

            /**
             * @brief Look up a value in the table based on indices. This is a constant-time operation. 
             */
            [[nodiscard]] double lookup(int q_index, int d_index) const override;

            /**
             * @brief Get the size of the table in the q-direction. 
             */
            [[nodiscard]] std::size_t size_q() const noexcept override;

            /**
             * @brief Get the size of the table in the d-direction. 
             */
            [[nodiscard]] std::size_t size_d() const noexcept override;

            /**
             * @brief Get an iterator to the beginning of the d-values for the given q-index.
             */
            [[nodiscard]] const constants::axes::d_type* begin(int q_index) const override;

            /**
             * @brief Get an iterator to the end of the d-values for the given q-index.
             */
            [[nodiscard]] const constants::axes::d_type* end(int q_index) const override;

            /**
             * @brief Get an iterator to the beginning of the d-values for the given q-index.
             */
            [[nodiscard]] constants::axes::d_type* begin(int q_index);

            /**
             * @brief Get an iterator to the end of the d-values for the given q-index.
             */
            [[nodiscard]] constants::axes::d_type* end(int q_index);

            /**
             * @brief Get the default table.
             */
            [[nodiscard]] static const VectorDebyeTable& get_default_table();

            /**
             * @brief Check if the two vectors are compatible with the default table. 
             *        Note that this check is only performed in debug mode.
             */
            static void check_default(const std::vector<double>& q, const std::vector<constants::axes::d_type>& d);

            /**
             * @brief Check if the vector is compatible with the default table. 
             *        Note that this check is only performed in debug mode.
             */
            static void check_default(const std::vector<constants::axes::d_type>& d);

        private: 
            /**
             * @brief Initialize this class for the given input. 
             * 
             * @param q The scattering vector to generate lookup values for. 
             * @param d The distance histogram to generate lookup values for. 
             */
            template<container_type T1, container_type T2>
            void initialize(const T1& q, const T2& d);

            /**
             * @brief Check if a DebyeTable is empty.
             */
            bool is_empty() const;
    };
}