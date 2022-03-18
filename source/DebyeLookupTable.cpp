#include <DebyeLookupTable.h>
#include <Utility.h>

using namespace table;

DebyeLookupTable::DebyeLookupTable() {}

DebyeLookupTable::DebyeLookupTable(const vector<double>& q, const vector<double>& d) {
    initialize(q, d);
}

void DebyeLookupTable::initialize(const vector<double>& q, const vector<double>& d) {
    // check if the default table can be used
    if (is_default(q, d)) {
        // check if the default table has already been instantiated                
        if (default_table.is_empty()) {
            double width = setting::axes::scattering_intensity_plot_binned_width;
            vector<double> _d(default_size/width, 0);
            for (unsigned int i = 1; i < d.size(); i++) {
                _d[i] = width*(i+0.5);
            }

            initialize(default_table, q, d);
        }

        // assign the lambda lookup function to a default table lookup
        lookup_function = [] (double q, double d) {return default_table.lookup(q, d);};
        index_lookup_function = [] (int i, int j) {return default_table.lookup_index(i, j);};
    } else {
        // assign the lambda lookup function to a custom table lookup
        initialize(table, q, d);
        lookup_function = [&] (double q, double d) {return table.lookup(q, d);};
        index_lookup_function = [&] (int i, int j) {return table.lookup_index(i, j);};
    }
}

/**
 * @brief Look up a value in the table based on @a q and @a d values. This is an amortized constant-time operation. 
 */
double DebyeLookupTable::lookup(double q, double d) const {
    return lookup_function(q, d);
}

/**
 * @brief Look up a value in the table based on indices. This is a constant-time operation. 
 */
double DebyeLookupTable::lookup(int row, int col) const {
    return index_lookup_function(row, col);
}

/**
 * @brief Check if this instance uses the default table. 
 */
bool DebyeLookupTable::uses_default_table() const {
    if (table.is_empty() && !default_table.is_empty()) {
        return true;
    }
    return false;
}

void DebyeLookupTable::initialize(LookupTable<double, double>& table, const vector<double>& q, const vector<double>& d) {
    table.initialize(q, d);
    for (const double Q : q) {
        for (const double D : d) {
            double qd = Q*D;
            double val = qd < tolerance ? 1 : sin(qd)/qd;
            table.assign(Q, D, val);
        }
    }
}

bool DebyeLookupTable::is_default(const vector<double>& q, const vector<double>& d) {
    // check q
    Axis& axis = setting::axes::scattering_intensity_plot_axis;
    double width = setting::axes::scattering_intensity_plot_binned_width;

    std::cout << "Checking if default table can be used.." << std::endl;
    if (q.size() != axis.bins) {return false;}
    if (q[0] != axis.min) {return false;}
    if (q[1] != axis.min + (axis.max-axis.min)/axis.bins) {return false;}
    if (q[2] != axis.min + 2*(axis.max-axis.min)/axis.bins) {return false;}

    // check d
    if (d[d.size()-1] > default_size) {return false;} // check if too large for default table
    if (!approx(d[2]-d[1], width)) {return false;} // check first width (d[1]-d[0] may be different from the default width)
    if (!approx(d[3]-d[2], width)) {return false;} // check second width

    return true;
}