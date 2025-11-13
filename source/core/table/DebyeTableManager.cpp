// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <table/DebyeTableManager.h>
#include <table/ArrayDebyeTable.h>

#include <stdexcept>
#include <numeric>
#include <cassert>

using namespace ausaxs;
using namespace ausaxs::table;

DebyeTableManager::DebyeTableManager() = default;
DebyeTableManager::DebyeTableManager(const DebyeTableManager& table) {*this = table;}
DebyeTableManager::DebyeTableManager(DebyeTableManager&&) noexcept = default; 
DebyeTableManager& DebyeTableManager::operator=(DebyeTableManager&&) noexcept = default;
DebyeTableManager& DebyeTableManager::operator=(const DebyeTableManager& table) {
    q = table.q;
    d = table.d;
    use_custom_table = table.use_custom_table;
    if (table.custom_sinc_table) {
        custom_sinc_table = std::make_unique<table::VectorDebyeTable>(d.axis, q.axis);
    }
    return *this;
}

observer_ptr<const table::DebyeTable> DebyeTableManager::get_sinc_table() const {
    if (use_custom_table) {
        if (recalculate) {
            if (q.defaulted && d.defaulted) {
                throw std::runtime_error("DebyeTableManager::get_sinc_table(): both q-axis and d-axis are defaulted, but custom table requested.");
            } else if (q.defaulted) {
                custom_sinc_table = std::make_unique<table::VectorDebyeTable>(d.axis, constants::axes::q_vals);
            } else if (d.defaulted) {
                custom_sinc_table = std::make_unique<table::VectorDebyeTable>(constants::axes::d_vals, q.axis);
            } else {
                custom_sinc_table = std::make_unique<table::VectorDebyeTable>(d.axis, q.axis);
            }
            recalculate = false;
        }
        return custom_sinc_table.get();
    } else {
        return &ArrayDebyeTable::get_default_table();
    }
}

void DebyeTableManager::reset_to_default() {
    q.defaulted = true;
    d.defaulted = true;
    use_custom_table = false;
}

bool appears_identical(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {return false;}
    if (a.size() > 10 && b.size() > 10) {
        for (size_t i = 0; i < 5; ++i) {
            if (a[i] != b[i]) {return false;}
        }
        for (size_t i = a.size()-5; i < a.size(); ++i) {
            if (a[i] != b[i]) {return false;}
        }
    }

    // extra assert in debug mode. if the above is not sufficient, this should catch it during testing
    assert(
        std::abs(std::accumulate(a.begin(), a.end(), 0.0) - std::accumulate(b.begin(), b.end(), 0.0)) < 1e-9 
        && "appears_identical: Sums do not match"
    );
    return true;
}

template<typename T, typename>
void DebyeTableManager::set_q_axis(T&& q_axis) {
    if (appears_identical(q.axis, q_axis)) {return;} // no change
    q.axis = std::forward<T>(q_axis);
    q.defaulted = false;
    use_custom_table = true;
    recalculate = true;
}

template<typename T, typename>
void DebyeTableManager::set_d_axis(T&& d_axis) {
    if (appears_identical(d.axis, d_axis)) {return;} // no change
    d.axis = std::forward<T>(d_axis);
    d.defaulted = false;
    use_custom_table = true;
    recalculate = true;
}

template void DebyeTableManager::set_q_axis(std::vector<double>&& q_axis);
template void DebyeTableManager::set_q_axis(const std::vector<double>& q_axis);
template void DebyeTableManager::set_q_axis(std::vector<double>& q_axis);
template void DebyeTableManager::set_d_axis(std::vector<double>&& d_axis);
template void DebyeTableManager::set_d_axis(const std::vector<double>& d_axis);
template void DebyeTableManager::set_d_axis(std::vector<double>& d_axis);