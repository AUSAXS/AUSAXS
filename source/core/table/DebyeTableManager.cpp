// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <table/DebyeTableManager.h>

#include <table/ArrayDebyeTable.h>
#include <utility/Exceptions.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::table;

DebyeTableManager::DebyeTableManager() = default;
DebyeTableManager::DebyeTableManager(const DebyeTableManager& table) {*this = table;}
DebyeTableManager::DebyeTableManager(DebyeTableManager&&) noexcept = default; 
DebyeTableManager& DebyeTableManager::operator=(DebyeTableManager&&) noexcept = default;
DebyeTableManager& DebyeTableManager::operator=(const DebyeTableManager& table) {
    if (this == &table) {return *this;}
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
                throw except::runtime_error("DebyeTableManager::get_sinc_table(): both q-axis and d-axis are defaulted, but custom table requested.");
            }
            
            if (q.defaulted) {
                custom_sinc_table = std::make_unique<table::VectorDebyeTable>(d.axis, constants::axes::q_vals);
            } else if (d.defaulted) {
                custom_sinc_table = std::make_unique<table::VectorDebyeTable>(constants::axes::d_vals, q.axis);
            } else {
                custom_sinc_table = std::make_unique<table::VectorDebyeTable>(d.axis, q.axis);
            }
            recalculate = false;
        }
        return custom_sinc_table.get();
    }
    return &ArrayDebyeTable::get_default_table();
}

void DebyeTableManager::reset_to_default() {
    q.defaulted = true;
    d.defaulted = true;
    use_custom_table = false;
}

template<typename T>
void DebyeTableManager::set_q_axis(T&& q_axis) requires (std::disjunction_v<
    std::is_rvalue_reference<T&&>,
    std::is_same<T, const std::vector<double>&>,
    std::is_same<T, std::vector<double>&>
>) {
    if (q_axis.size() <= q.axis.size() && std::equal(q_axis.begin(), q_axis.end(), q.axis.begin())) {return;}
    q.axis = std::forward<T>(q_axis);
    q.defaulted = false;
    use_custom_table = true;
    recalculate = true;
}

template<typename T>
void DebyeTableManager::set_d_axis(T&& d_axis) requires (std::disjunction_v<
    std::is_rvalue_reference<T&&>,
    std::is_same<T, const std::vector<double>&>,
    std::is_same<T, std::vector<double>&>
>) {
    if (d_axis.size() == d.axis.size() && std::equal(d_axis.begin(), d_axis.end(), d.axis.begin())) {return;}
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