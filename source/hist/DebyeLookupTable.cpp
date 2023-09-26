#include <hist/DebyeLookupTable.h>
#include <utility/Console.h>
#include <utility/Utility.h>
#include <utility/Axis.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>

#include <cmath>

using namespace table;

DebyeLookupTable::DebyeLookupTable(const std::vector<double>& q, const std::vector<double>& d) {
    initialize(q, d);
}

void DebyeLookupTable::initialize(LookupTable<double, double>& table, const std::vector<double>& q, const std::vector<double>& d) {
    table.initialize(q, d);
    for (unsigned int i = 0; i < q.size(); ++i) {
        for (unsigned int j = 0; j < d.size(); ++j) {
            double qd = q[i]*d[j];  
            double val;
            if (qd < tolerance) [[unlikely]] {
                double qd2 = qd*qd;
                val = 1 - qd2/6 + qd2*qd2/120;
            } else {
                val = std::sin(qd)/qd;
            }
            table.assign_index(i, j, val);
        }
    }
}

void DebyeLookupTable::initialize(const std::vector<double>& q, const std::vector<double>& d) {
    // check if the default table can be used
    if (is_default(q, d)) {
        // check if the default table has already been instantiated 
        if (default_table.is_empty()) {
            std::vector<double> _d = Axis(0, settings::axes::max_distance, settings::axes::max_distance/settings::axes::distance_bin_width).as_vector(0.5);
            _d[0] = 0; // fix the first bin to 0 since it primarily contains self-correlation terms
            initialize(default_table, q, _d); // note we pass _d and not d
        }

        // assign the lambda lookup function to a default table lookup
        lookup_function = [] (double q, double d) {return default_table.lookup(q, d);};
        index_lookup_function = [] (unsigned int i, unsigned int j) {return default_table.lookup_index(i, j);};
    } else {
        if (settings::general::verbose) {console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables. ");}

        // assign the lambda lookup function to a custom table lookup
        initialize(table, q, d);
        lookup_function = [this] (double q, double d) {return table.lookup(q, d);};
        index_lookup_function = [this] (unsigned int i, unsigned int j) {return table.lookup_index(i, j);};
    }
}

DebyeLookupTable DebyeLookupTable::get_default_table() {
    std::vector<double> _d = Axis(0, settings::axes::max_distance, settings::axes::max_distance/settings::axes::distance_bin_width).as_vector(0.5);
    _d[0] = 0;
    std::vector<double> _q = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
    return DebyeLookupTable(_q, _d);
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
double DebyeLookupTable::lookup(unsigned int q_index, unsigned int d_index) const {
    return index_lookup_function(q_index, d_index);
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

#define NOT_DEFAULT_MSG false
#if NOT_DEFAULT_MSG
    #include <iostream>
#endif
bool DebyeLookupTable::is_default(const std::vector<double>& q, const std::vector<double>& d) {
    // check q
    Axis axis = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);
    double width = settings::axes::distance_bin_width;

    if (q.size() != axis.bins) [[unlikely]] {
        #if NOT_DEFAULT_MSG
            std::cout << "q.size() != axis.bins" << std::endl;
        #endif
        return false;
    }

    if (q[0] != axis.min) [[unlikely]] {
        #if NOT_DEFAULT_MSG
            std::cout << "q[0] != axis.min" << std::endl;
        #endif
        return false;
    }

    if (q[1] != axis.min + (axis.max-axis.min)/axis.bins) [[unlikely]] {
        #if NOT_DEFAULT_MSG
            std::cout << "q[1] != axis.min + (axis.max-axis.min)/axis.bins" << std::endl;
        #endif
        return false;
    }

    if (q[2] != axis.min + 2*(axis.max-axis.min)/axis.bins) [[unlikely]] {
        #if NOT_DEFAULT_MSG
            std::cout << "q[2] != axis.min + 2*(axis.max-axis.min)/axis.bins" << std::endl;
        #endif        
        return false;
    }

    // check empty
    if (d.empty()) [[unlikely]] {
        #if NOT_DEFAULT_MSG
            std::cout << "d.empty()" << std::endl;
        #endif
        return false;
    }

    // check if too large for default table
    if (d.back() > settings::axes::max_distance) [[unlikely]] {
        #if NOT_DEFAULT_MSG
            std::cout << "d.back() > default_size" << std::endl;
        #endif
        return false;
    }
    
    // check first width (d[1]-d[0] may be different from the default width)
    if (!utility::approx(d[2]-d[1], width)) [[unlikely]] {
        #if NOT_DEFAULT_MSG
            std::cout << "!utility::approx(d[2]-d[1], width)" << std::endl;
        #endif
        return false;
    }
    
    // check second width
    if (!utility::approx(d[3]-d[2], width)) [[unlikely]] {
        #if NOT_DEFAULT_MSG
            std::cout << "!utility::approx(d[3]-d[2], width)" << std::endl;
        #endif
        return false;
    }

    return true;
}

void DebyeLookupTable::reset_default_table() {
    default_table = LookupTable<double, double>();
}

unsigned int DebyeLookupTable::size_q() const {
    return table.N;
}

unsigned int DebyeLookupTable::size_d() const {
    return table.M;
}