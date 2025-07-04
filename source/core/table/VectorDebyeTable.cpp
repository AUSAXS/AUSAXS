// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <table/VectorDebyeTable.h>
#include <utility/Console.h>
#include <utility/Utility.h>
#include <utility/Axis.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <constants/Constants.h>
#include <math/ConstexprMath.h>

#include <cmath>

using namespace ausaxs;
using namespace ausaxs::table;

VectorDebyeTable::VectorDebyeTable() = default;

VectorDebyeTable::VectorDebyeTable(const std::vector<constants::axes::d_type>& d) : Table(constants::axes::q_axis.bins, d.size()) {
    initialize(constants::axes::q_vals, d);
}

VectorDebyeTable::VectorDebyeTable(const std::array<constants::axes::d_type, constants::axes::d_axis.bins>& d) : Table(constants::axes::q_axis.bins, d.size()) {
    initialize(constants::axes::q_vals, d);
}

VectorDebyeTable::VectorDebyeTable(const std::vector<constants::axes::d_type>& d, const std::vector<double>& q) : Table(q.size(), d.size()) {
    initialize(q, d);
}

template<container_type T1, container_type T2>
void VectorDebyeTable::initialize(const T1& q, const T2& d) {
    constexpr double tolerance = 1e-3;  // The minimum x-value where sin(x)/x is replaced by its Taylor-series.
    constexpr double inv_6 = 1./6;      // 1/6
    constexpr double inv_120 = 1./120;  // 1/120

    for (int i = 0; i < static_cast<int>(N); ++i) {
        for (int j = 0; j < static_cast<int>(M); ++j) {
            double qd = q[i]*d[j];
            if (qd < tolerance) {
                double qd2 = qd*qd;
                index(i, j) = 1 - qd2*inv_6 + qd2*qd2*inv_120;
            } else {
                // index(i, j) = constexpr_math::fast::sin(qd)/qd;
                index(i, j) = std::sin(qd)/qd;
            }
        }
    }
}

bool VectorDebyeTable::is_empty() const {return N == 0 || M == 0;}

double VectorDebyeTable::lookup(int q_index, int d_index) const {
    return index(q_index, d_index);
}

const constants::axes::d_type* VectorDebyeTable::begin(int q_index) const {
    return data.data() + q_index*M;
}

const constants::axes::d_type* VectorDebyeTable::end(int q_index) const {
    return data.data() + (q_index+1)*M;
}

constants::axes::d_type* VectorDebyeTable::begin(int q_index) {
    return data.data() + q_index*M;
}

constants::axes::d_type* VectorDebyeTable::end(int q_index) {
    return data.data() + (q_index+1)*M;
}

std::size_t VectorDebyeTable::size_q() const noexcept {
    return N;
}

std::size_t VectorDebyeTable::size_d() const noexcept {
    return M;
}

const VectorDebyeTable& VectorDebyeTable::get_default_table() {
    static VectorDebyeTable default_table(constants::axes::d_vals);
    return default_table;
}

#if DEBUG 
    #include <iostream>
    #include <utility/Console.h>
    #include <utility/Utility.h>
    #include <settings/HistogramSettings.h>
    #include <settings/GeneralSettings.h>
#endif
void VectorDebyeTable::check_default(const std::vector<double>& q, const std::vector<constants::axes::d_type>& d) {
    #if DEBUG 
        if (!settings::general::warnings) {return;}
        const Axis& axis = constants::axes::q_axis;

        auto qvals = axis.as_vector();
        unsigned int i = 0;
        for (; i < axis.bins; ++i) {
            if (utility::approx(q.front(), qvals[i])) {break;}
        }
        if (i == axis.bins) [[unlikely]] {
            console::print_warning("Warning in VectorDebyeTable::initialize: Incompatible with default tables.");
            std::cout << "\tReason: q[0] does not match any index of default q-array" << std::endl;

        }
        if (q[0] != qvals[i]) [[unlikely]] {
            console::print_warning("Warning in VectorDebyeTable::initialize: Incompatible with default tables.");
            std::cout << "\tReason: q[0] != axis.min" << std::endl;
        }

        if (q[1] != qvals[i+1]) [[unlikely]] {
            console::print_warning("Warning in VectorDebyeTable::initialize: Incompatible with default tables.");
            std::cout << "\tReason: q[1] != axis.min + (axis.max-axis.min)/axis.bins" << std::endl;
        }

        if (q[2] != qvals[i+2]) [[unlikely]] {
            console::print_warning("Warning in VectorDebyeTable::initialize: Incompatible with default tables.");
            std::cout << "\tReason: q[2] != axis.min + 2*(axis.max-axis.min)/axis.bins" << std::endl;
        }

        check_default(d);
    #endif
}

void VectorDebyeTable::check_default(const std::vector<constants::axes::d_type>& d) {
    #ifdef DEBUG
        // check empty
        if (d.empty()) [[unlikely]] {
            console::print_warning("Warning in VectorDebyeTable::initialize: Incompatible with default tables.");
            std::cout << "\tReason: d.empty()" << std::endl;
        }

        // check if too large for default table
        if (d.back() > constants::axes::d_axis.max) [[unlikely]] {
            console::print_warning("Warning in VectorDebyeTable::initialize: Incompatible with default tables.");
            std::cout << "\tReason: d.back() > default_size" << std::endl;
        }
        
        // check first width (d[1]-d[0] may be different from the default width)
        if (!utility::approx(d[2]-d[1], constants::axes::d_axis.width())) [[unlikely]] {
            console::print_warning("Warning in VectorDebyeTable::initialize: Incompatible with default tables.");
            std::cout << "\tReason: !utility::approx(d[2]-d[1], width)" << std::endl;
        }
        
        // check second width
        if (!utility::approx(d[3]-d[2], constants::axes::d_axis.width())) [[unlikely]] {
            console::print_warning("Warning in VectorDebyeTable::initialize: Incompatible with default tables.");
            std::cout << "\tReason: !utility::approx(d[3]-d[2], width)" << std::endl;
        }
    #endif
}