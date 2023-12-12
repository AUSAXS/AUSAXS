#include <table/VectorDebyeTable.h>
#include <utility/Console.h>
#include <utility/Utility.h>
#include <utility/Axis.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <constants/Constants.h>
#include <math/ConstexprMath.h>

#include <cmath>

using namespace table;

VectorDebyeTable::VectorDebyeTable() = default;

VectorDebyeTable::VectorDebyeTable(const std::vector<constants::axes::d_type>& d) : Table(constants::axes::q_axis.bins, d.size()) {
    initialize(constants::axes::q_vals, d);
}

VectorDebyeTable::VectorDebyeTable(const std::array<constants::axes::d_type, constants::axes::d_axis.bins>& d) : Table(constants::axes::q_axis.bins, d.size()) {
    initialize(constants::axes::q_vals, d);
}

template<container_type T1, container_type T2>
void VectorDebyeTable::initialize(const T1& q, const T2& d) {
    constexpr double tolerance = 1e-3;  // The minimum x-value where sin(x)/x is replaced by its Taylor-series.
    constexpr double inv_6 = 1./6;      // 1/6
    constexpr double inv_120 = 1./120;  // 1/120

    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < M; ++j) {
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

/**
 * @brief Look up a value in the table based on indices. This is a constant-time operation. 
 */
double VectorDebyeTable::lookup(unsigned int q_index, unsigned int d_index) const {
    return index(q_index, d_index);
}

const constants::axes::d_type* VectorDebyeTable::begin(unsigned int q_index) const {
    return data.data() + q_index*M;
}

const constants::axes::d_type* VectorDebyeTable::end(unsigned int q_index) const {
    return data.data() + (q_index+1)*M;
}

constants::axes::d_type* VectorDebyeTable::begin(unsigned int q_index) {
    return data.data() + q_index*M;
}

constants::axes::d_type* VectorDebyeTable::end(unsigned int q_index) {
    return data.data() + (q_index+1)*M;
}

unsigned int VectorDebyeTable::size_q() const noexcept {
    return N;
}

unsigned int VectorDebyeTable::size_d() const noexcept {
    return M;
}