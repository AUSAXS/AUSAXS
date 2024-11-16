/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#ifdef CONSTEXPR_LOOKUP_TABLE
    #include <table/ArrayDebyeTable.h>

    using namespace ausaxs;
    using namespace ausaxs::table;

    #if DEBUG 
        #include <iostream>
        #include <utility/Console.h>
        #include <utility/Utility.h>
        #include <settings/HistogramSettings.h>
        #include <settings/GeneralSettings.h>
    #endif
    void ArrayDebyeTable::check_default(const std::vector<double>& q, const std::vector<double>& d) {
        #if DEBUG 
            if (!settings::general::warnings) {return;}
            const Axis& axis = constants::axes::q_axis;

            auto qvals = axis.as_vector();
            unsigned int i = 0;
            for (; i < axis.bins; ++i) {
                if (utility::approx(q.front(), qvals[i])) {break;}
            }
            if (i == axis.bins) [[unlikely]] {
                console::print_warning("Warning in ArrayDebyeTable::initialize: Incompatible with default tables.");
                std::cout << "\tReason: q[0] does not match any index of default q-array" << std::endl;

            }
            if (q[0] != qvals[i]) [[unlikely]] {
                console::print_warning("Warning in ArrayDebyeTable::initialize: Incompatible with default tables.");
                std::cout << "\tReason: q[0] != axis.min" << std::endl;
            }

            if (q[1] != qvals[i+1]) [[unlikely]] {
                console::print_warning("Warning in ArrayDebyeTable::initialize: Incompatible with default tables.");
                std::cout << "\tReason: q[1] != axis.min + (axis.max-axis.min)/axis.bins" << std::endl;
            }

            if (q[2] != qvals[i+2]) [[unlikely]] {
                console::print_warning("Warning in ArrayDebyeTable::initialize: Incompatible with default tables.");
                std::cout << "\tReason: q[2] != axis.min + 2*(axis.max-axis.min)/axis.bins" << std::endl;
            }

            check_default(d);
        #endif
    }

    void ArrayDebyeTable::check_default(const std::vector<constants::axes::d_type>& d) {
        #ifdef DEBUG
            // check empty
            if (d.empty()) [[unlikely]] {
                console::print_warning("Warning in ArrayDebyeTable::initialize: Incompatible with default tables.");
                std::cout << "\tReason: d.empty()" << std::endl;
            }

            // check if too large for default table
            if (d.back() > constants::axes::d_axis.max) [[unlikely]] {
                console::print_warning("Warning in ArrayDebyeTable::initialize: Incompatible with default tables.");
                std::cout << "\tReason: d.back() > default_size" << std::endl;
            }
            
            // check first width (d[1]-d[0] may be different from the default width)
            if (!utility::approx(d[2]-d[1], constants::axes::d_axis.width())) [[unlikely]] {
                console::print_warning("Warning in ArrayDebyeTable::initialize: Incompatible with default tables.");
                std::cout << "\tReason: !utility::approx(d[2]-d[1], width)" << std::endl;
            }
            
            // check second width
            if (!utility::approx(d[3]-d[2], constants::axes::d_axis.width())) [[unlikely]] {
                console::print_warning("Warning in ArrayDebyeTable::initialize: Incompatible with default tables.");
                std::cout << "\tReason: !utility::approx(d[3]-d[2], width)" << std::endl;
            }
        #endif
    }

    #pragma message("Precompiling sinc lookup table. This may take a couple of minutes...")
    
    inline constexpr ArrayDebyeTable default_table;
    const ArrayDebyeTable& ArrayDebyeTable::get_default_table() {
        return default_table;
    }
#endif