// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/distribution/detail/WeightedEntry.h>

#include <iostream>

using namespace ausaxs;
using namespace ausaxs::hist::detail;

std::ostream& hist::detail::operator<<(std::ostream& os, const WeightedEntry& entry) {
    os << "WeightedEntry(" << entry.value << ", " << entry.count << ", " << entry.bin_center << ")";
    return os;
}