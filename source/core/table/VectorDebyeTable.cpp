// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <table/VectorDebyeTable.h>

#include <constants/ConstantsAxes.h>
#include <math/ConstexprMath.h>
#include <utility/MultiThreading.h>

using namespace ausaxs;
using namespace ausaxs::table;

VectorDebyeTable::VectorDebyeTable() = default;

template<container_type T1, container_type T2>
VectorDebyeTable::VectorDebyeTable(const T1& d, const T2& q) : Table(q.size(), d.size()) {initialize(q, d);}

template<container_type T1, container_type T2>
void VectorDebyeTable::initialize(const T1& q, const T2& d) {
    auto* pool = utility::multi_threading::get_global_pool();
    pool->submit_blocks(0, N, // iterate through all q values
        [this, &q, &d] (int start, int end) {
            for (int i = start; i < end; ++i) {
                for (int j = 0; j < M; ++j) {
                    index(i, j) = math::sinc(q[i]*d[j]);
                }
            }
        }
    ).wait();
}

template VectorDebyeTable::VectorDebyeTable(const std::vector<double>&, const std::vector<double>&);
template VectorDebyeTable::VectorDebyeTable(const std::vector<double>&, const std::array<double, constants::axes::q_axis.bins>&);
template VectorDebyeTable::VectorDebyeTable(const std::array<double, constants::axes::d_axis.bins>&, const std::vector<double>&);
template VectorDebyeTable::VectorDebyeTable(const std::array<double, constants::axes::d_axis.bins>&, const std::array<double, constants::axes::q_axis.bins>&);

bool VectorDebyeTable::is_empty() const {return N == 0 || M == 0;}

double VectorDebyeTable::lookup(int q_index, int d_index) const {
    return index(q_index, d_index);
}

const double* VectorDebyeTable::begin(int q_index) const {
    return data.data() + q_index*M;
}

const double* VectorDebyeTable::end(int q_index) const {
    return data.data() + (q_index+1)*M;
}

double* VectorDebyeTable::begin(int q_index) {
    return data.data() + q_index*M;
}

double* VectorDebyeTable::end(int q_index) {
    return data.data() + (q_index+1)*M;
}

std::size_t VectorDebyeTable::size_q() const noexcept {
    return N;
}

std::size_t VectorDebyeTable::size_d() const noexcept {
    return M;
}

const VectorDebyeTable& VectorDebyeTable::get_default_table() {
    static VectorDebyeTable default_table(constants::axes::d_vals, constants::axes::q_vals);
    return default_table;
}

