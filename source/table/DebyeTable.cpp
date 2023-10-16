#include <table/DebyeTable.h>
#include <utility/Console.h>
#include <utility/Utility.h>
#include <utility/Axis.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>

#include <cmath>

using namespace table;

DebyeTable::DebyeTable() = default;

DebyeTable::DebyeTable(const std::vector<double>& q, const std::vector<double>& d) : Table(q.size(), d.size()) {
    initialize(q, d);
}

// fast cosine from https://stackoverflow.com/a/28050328
// error margin is 0.00109 in the range [-pi, pi]
[[maybe_unused]] inline double fast_cos(double x) noexcept {
    constexpr double tp = 1./(2.*M_PI);
    x *= tp;
    x -= 0.25 + std::floor(x + 0.25);
    x *= 16.*(std::abs(x) - 0.5);
    x += .225*x*(std::abs(x) - 1.);
    return x;
}

// conversion of above expression but for sine
inline double fast_sin(double x) {
    constexpr double tp = 1./(2.*M_PI);
    x *= tp;
    x -= 0.5 + std::floor(x);
    x *= 16.*(std::abs(x) - 0.5);
    x += .225*x*(std::abs(x) - 1.);
    return x;
}

void DebyeTable::initialize(const std::vector<double>& q, const std::vector<double>& d) {
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
                index(i, j) = fast_sin(qd)/qd;
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
        Axis axis = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);
        double width = settings::axes::distance_bin_width;

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
        if (d.back() > settings::axes::max_distance) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: d.back() > default_size" << std::endl;
        }
        
        // check first width (d[1]-d[0] may be different from the default width)
        if (!utility::approx(d[2]-d[1], width)) [[unlikely]] {
            console::print_warning("Warning in DebyeLookupTable::initialize: Not using default tables.");
            std::cout << "\tReason: !utility::approx(d[2]-d[1], width)" << std::endl;
        }
        
        // check second width
        if (!utility::approx(d[3]-d[2], width)) [[unlikely]] {
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
        std::vector<double> _d = Axis(0, settings::axes::max_distance, settings::axes::max_distance/settings::axes::distance_bin_width).as_vector(0.5);
        _d[0] = 0; // fix the first bin to 0 since it primarily contains self-correlation terms
        std::vector<double> _q = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
        default_table = DebyeTable(_q, _d);
    }
    return default_table;
}

void DebyeTable::reset_default_table() {
    default_table = DebyeTable();
}