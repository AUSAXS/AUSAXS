#include <table/DebyeTable.h>
#include <utility/Console.h>
#include <utility/Utility.h>
#include <utility/Axis.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <constants/Constants.h>
#include <math/Trigonometry.h>

#include <cmath>

using namespace table;

DebyeTable::DebyeTable() = default;

template<container_type T1, container_type T2>
DebyeTable::DebyeTable(const T1& q, const T2& d) : Table(q.size(), d.size()) {
    initialize(q, d);
}

template<container_type T1, container_type T2>
void DebyeTable::initialize(const T1& q, const T2& d) {
    constexpr double tolerance = 1;     // The minimum x-value where sin(x)/x is replaced by its Taylor-series.
    constexpr double inv_6 = 1./6;      // 1/6
    constexpr double inv_120 = 1./120;  // 1/120

    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < M; ++j) {
            double qd = q[i]*d[j];
            if (qd < tolerance) {
                double qd2 = qd*qd;
                index(i, j) = 1 - qd2*inv_6 + qd2*qd2*inv_120;
            } else {
                index(i, j) = math::fast::sin(qd)/qd;
            }
        }
    }
}

bool DebyeTable::is_empty() const {return N == 0 || M == 0;}

/**
 * @brief Look up a value in the table based on indices. This is a constant-time operation. 
 */
double DebyeTable::lookup(unsigned int q_index, unsigned int d_index) const {
    return index(q_index, d_index);
}

const std::vector<double>::const_iterator DebyeTable::begin(unsigned int q_index) const {
    return Table::begin(q_index);
}

const std::vector<double>::const_iterator DebyeTable::end(unsigned int q_index) const {
    return Table::end(q_index);
}

std::vector<double>::iterator DebyeTable::begin(unsigned int q_index) {
    return Table::begin(q_index);
}

std::vector<double>::iterator DebyeTable::end(unsigned int q_index) {
    return Table::end(q_index);
}

#if DEBUG 
    #include <iostream>
#endif
void DebyeTable::check_default([[maybe_unused]] const std::vector<double>& q, [[maybe_unused]] const std::vector<double>& d) {
    #if DEBUG
        Axis axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);

        if (q.size() != axis.bins) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: q.size() != axis.bins" << std::endl;
        }

        if (q[0] != axis.min) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: q[0] != axis.min" << std::endl;
        }

        if (q[1] != axis.min + (axis.max-axis.min)/axis.bins) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: q[1] != axis.min + (axis.max-axis.min)/axis.bins" << std::endl;
        }

        if (q[2] != axis.min + 2*(axis.max-axis.min)/axis.bins) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: q[2] != axis.min + 2*(axis.max-axis.min)/axis.bins" << std::endl;
        }

        // check empty
        if (d.empty()) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: d.empty()" << std::endl;
        }

        // check if too large for default table
        if (d.back() > constants::axes::d_axis.max) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: d.back() > default_size" << std::endl;
        }
        
        // check first width (d[1]-d[0] may be different from the default width)
        if (!utility::approx(d[2]-d[1], constants::axes::d_axis.width())) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: !utility::approx(d[2]-d[1], width)" << std::endl;
        }
        
        // check second width
        if (!utility::approx(d[3]-d[2], constants::axes::d_axis.width())) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: !utility::approx(d[3]-d[2], width)" << std::endl;
        }
    #endif
}

unsigned int DebyeTable::size_q() const {
    return N;
}

unsigned int DebyeTable::size_d() const {
    return M;
}

DebyeTable DebyeTable::default_table = DebyeTable();
const DebyeTable& DebyeTable::get_default_table() {
    if (default_table.is_empty()) {
        auto& d = constants::axes::d_vals;
        auto& q = constants::axes::q_vals;
        default_table = DebyeTable(q, d);
    }
    return default_table;
}

void DebyeTable::reset_default_table() {
    default_table = DebyeTable();
}