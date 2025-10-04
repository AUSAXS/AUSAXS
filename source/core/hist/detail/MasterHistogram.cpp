// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/MasterHistogram.h>

using namespace ausaxs;
using namespace ausaxs::hist::detail;

template<bool use_weighted_distribution> 
MasterHistogram<use_weighted_distribution>::MasterHistogram() = default;

template<bool use_weighted_distribution> 
MasterHistogram<use_weighted_distribution>::MasterHistogram(std::vector<double>&& p_base, const Axis& axis) : GenericDistribution1D_t(p_base), base(std::move(p_base)), axis(axis) {}

template<bool use_weighted_distribution> 
MasterHistogram<use_weighted_distribution>::MasterHistogram(const std::vector<double>& p_base, const Axis& axis) : GenericDistribution1D_t(p_base), base(p_base), axis(axis) {}

template<>
MasterHistogram<true>& MasterHistogram<true>::operator+=(const GenericDistribution1D_t& rhs) {
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<>());
    return *this; 
}

template<bool use_weighted_distribution>
MasterHistogram<use_weighted_distribution>& MasterHistogram<use_weighted_distribution>::operator+=(const GenericDistribution1D_t& rhs) {
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<>());
    return *this; 
}

template<bool use_weighted_distribution> 
MasterHistogram<use_weighted_distribution>& MasterHistogram<use_weighted_distribution>::operator-=(const GenericDistribution1D_t& rhs) {
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::minus<>());
    return *this;
}

template class hist::detail::MasterHistogram<true>;
template class hist::detail::MasterHistogram<false>;