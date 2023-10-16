#pragma once

#include <table/LookupTable.h>

#include <functional>

namespace table {
    /**
     * @brief Specialized class for storing the sin(x)/x lookup table as required in the Debye equation. 
     *        The primary feature of this class is that it checks if a default table can be used for a given instance, and if so, 
     *        calculates and uses only a single table for multiple instances. The default table is static and can thus be used for
     *        later instances even should the first one be destroyed.
     * 
     * If the default table cannot be used, this class works exactly like the standard \class LookupTable. 
     */
    class [[deprecated("use DebyeTable instead")]] DebyeLookupTable {
        public: 
            /**
             * @brief Constructor. 
             * 
             * @param q The scattering vector to generate lookup values for. 
             * @param d The distance histogram to generate lookup values for. 
             */
            DebyeLookupTable(const std::vector<double>& q, const std::vector<double>& d);

            /**
             * @brief Look up a value in the table based on @a q and @a d values. This is an amortized constant-time operation. 
             */
            [[nodiscard]] double lookup(double q, double d) const;

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
             * @brief Check if this instance uses the default table. 
             */
            [[nodiscard]] bool uses_default_table() const;

            /**
             * @brief Reset the default table. 
             */
            static void reset_default_table();

            /**
             * @brief Get the default table. 
             */
            [[nodiscard]] static DebyeLookupTable get_default_table();

        private: 
            std::function<double(double, double)> lookup_function;                      // The lookup function when providing q and d values. 
            std::function<double(unsigned int, unsigned int)> index_lookup_function;    // The lookup function when providing indices in the table. 
            LookupTable<double, double> table;                                          // The specialized table for this specific instance. 
            inline static LookupTable<double, double> default_table;                    // The shared default table. 
            inline static double tolerance = 1e-6;                                      // The minimum x-value where sin(x)/x is replaced by its Taylor-series.

            /**
             * @brief Initialize this class for the given input. 
             * 
             * @param q The scattering vector to generate lookup values for. 
             * @param d The distance histogram to generate lookup values for. 
             */
            void initialize(const std::vector<double>& q, const std::vector<double>& d);

            /**
             * @brief Initialize a given table based on the input vectors.  
             * 
             * @param table The table to be initialized. Should either be the table or default_table members of this class. 
             * @param q The scattering vector to generate lookup values for. 
             * @param d The distance histogram to generate lookup values for. 
             */
            static void initialize(LookupTable<double, double>& table, const std::vector<double>& q, const std::vector<double>& d);

            /**
             * @brief Determine if the two arguments are based on the default settings. 
             * 
             * @param q The scattering vector.
             * @param d The histogram bins.
             */
            [[nodiscard]] static bool is_default(const std::vector<double>& q, const std::vector<double>& d);
    };
}