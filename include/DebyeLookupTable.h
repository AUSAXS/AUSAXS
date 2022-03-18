#include <settings.h>
#include <Table.h>

namespace table {
    class DebyeLookupTable {
        public: 
            DebyeLookupTable();

            DebyeLookupTable(const vector<double>& q, const vector<double>& d);

            void initialize(const vector<double>& q, const vector<double>& d);

            /**
             * @brief Look up a value in the table based on @a q and @a d values. This is an amortized constant-time operation. 
             */
            double lookup(double q, double d) const;

            /**
             * @brief Look up a value in the table based on indices. This is a constant-time operation. 
             */
            double lookup(int row, int col) const;

            /**
             * @brief Check if this instance uses the default table. 
             */
            bool uses_default_table() const;

        private: 
            std::function<double(double, double)> lookup_function;
            std::function<double(int, int)> index_lookup_function;
            LookupTable<double, double> table;
            inline static LookupTable<double, double> default_table;
            inline static int default_size = 500; 
            inline static double tolerance = 1e-9;

            static void initialize(LookupTable<double, double>& table, const vector<double>& q, const vector<double>& d);

            /**
             * @brief Determine if the two arguments are based on the default settings. 
             * 
             * @param q The scattering vector
             * @param d The histogram bins
             */
            static bool is_default(const vector<double>& q, const vector<double>& d);
    };
}