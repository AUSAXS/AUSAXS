#include <vector>
#include <unordered_map>
#include <iostream>

using std::vector;

class Table {
    public: 
        Table() : N(0), M(0), data(0) {}

        Table(unsigned int N, unsigned int M) : N(N), M(M), data(N*M) {}

        const double& index(unsigned int i, unsigned int j) const {return data[M*i + j];}

        double& index(unsigned int i, unsigned int j) {return data[M*i + j];}

        unsigned int N, M;
        vector<double> data;
};

template<typename T, typename Q>
class LookupTable : public Table {
    public: 
        LookupTable() {}

        LookupTable(const vector<T>& rows, const vector<Q>& columns) {
            initialize(rows, columns);
        }

        void initialize(const vector<T>& rows, const vector<Q>& columns) {
            N = rows.size(); M = columns.size();
            data.resize(N*M);

            Tmap.reserve(N);
            for (unsigned int i = 0; i < N; i++) {
                Tmap[rows[i]] = i;
            }

            Qmap.reserve(M);
            for (unsigned int i = 0; i < M; i++) {
                Qmap[columns[i]] = i;
            }
        }

        bool is_empty() const {return data.size() == 0;}

        double lookup(const T row, const Q col) const {return index(Tmap[row], Qmap[col]);}

        void assign(const T row, const Q col, double val) {index(Tmap[row], Qmap[col]) = val;}

    private:
        std::unordered_map<T, unsigned int> Tmap;        
        std::unordered_map<Q, unsigned int> Qmap;        
};

#include <settings.h>
class DebyeLookupTable {
    public: 
        DebyeLookupTable() {}

        void initialize(const vector<double>& q, const vector<double>& d) {
            // check if the default table can be used
            if (is_default(q, d)) {

                // check if the default table has already been instantiated                
                if (default_table.is_empty()) {
                    double width = setting::axes::scattering_intensity_plot_binned_width;
                    vector<double> _d(default_size/width, 0);
                    for (int i = 1; i < d.size(); i++) {
                        _d[i] = width*(i+0.5);
                    }

                    default_table.initialize(q, _d);
                }

                // assign the lambda lookup function to a default table lookup
                lookup_function = [] (double q, double d) {return default_table.lookup(q, d);};
            } else {

                // assign the lambda lookup function to a custom table lookup
                table.initialize(q, d);
                lookup_function = [&] (double q, double d) {return table.lookup(q, d);};
            }
        }

        double lookup(double q, double d) const {
            return lookup_function(q, d);
        }

        /**
         * @brief Determine if the two arguments are based on the default settings. 
         * 
         * @param q The scattering vector
         * @param d The histogram bins
         */
        bool is_default(const vector<double>& q, const vector<double>& d) const {
            // check q
            if (q.size() != setting::axes::scattering_intensity_plot_axis.bins) {return false;}
            if (q[0] != setting::axes::scattering_intensity_plot_axis.min) {return false;}
            if (q[q.size()-1] != setting::axes::scattering_intensity_plot_axis.max) {return false;}

            // check d
            if (d[d.size()-1] > default_size) {return false;} // check if too large for default table
            if (d[2]-d[1] != setting::axes::scattering_intensity_plot_binned_width) {return false;} // check first width (d[1]-d[0] may be different from the default width)
            if (d[3]-d[2] != setting::axes::scattering_intensity_plot_binned_width) {return false;} // check second width

            return true;
        }

    private: 
        std::function<double(double, double)> lookup_function;
        LookupTable<double, double> table;
        static LookupTable<double, double> default_table;
        inline static int default_size = 500; 
};