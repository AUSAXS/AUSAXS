#pragma once

#include <table/Table.h>
#include <utility/Concepts.h>

namespace table {
    class DebyeTable : private Table {
        public:
            DebyeTable();

            /**
             * @brief Constructor. 
             * 
             * @param q The scattering vector to generate lookup values for. 
             * @param d The distance histogram to generate lookup values for. 
             */
            template<container_type T1, container_type T2>
            DebyeTable(const T1& q, const T2& d);

            /**
             * @brief Look up a value in the table based on indices. This is a constant-time operation. 
             */
            [[nodiscard]] double lookup(unsigned int q_index, unsigned int d_index) const;

            /**
             * @brief Get the size of the table in the q-direction. 
             */
            [[nodiscard]] unsigned int size_q() const;

            /**
             * @brief Get the size of the table in the d-direction. 
             */
            [[nodiscard]] unsigned int size_d() const;

            /**
             * @brief Get an iterator to the beginning of the d-values for the given q-index.
             */
            const std::vector<double>::const_iterator begin(unsigned int q_index) const;

            /**
             * @brief Get an iterator to the end of the d-values for the given q-index.
             */
            const std::vector<double>::const_iterator end(unsigned int q_index) const;

            /**
             * @brief Get an iterator to the beginning of the d-values for the given q-index.
             */
            std::vector<double>::iterator begin(unsigned int q_index);

            /**
             * @brief Get an iterator to the end of the d-values for the given q-index.
             */
            std::vector<double>::iterator end(unsigned int q_index);

            /**
             * @brief Reset the default table. 
             */
            static void reset_default_table();

            /**
             * @brief Get the default table. 
             */
            [[nodiscard]] static const DebyeTable& get_default_table();

            /**
             * @brief Check if the two vectors are compatible with the default table. 
             *        Note that this check is only performed in debug mode.
             */
            static void check_default(const std::vector<double>& q, const std::vector<double>& d);

        private: 
            static DebyeTable default_table;

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