/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/distribution/detail/WeightedEntry.h>

#include <iostream>

using namespace ausaxs;
using namespace ausaxs::hist::detail;

std::ostream& hist::detail::operator<<(std::ostream& os, const WeightedEntry& entry) {
    os << "WeightedEntry(" << entry.value << ", " << entry.count << ", " << entry.bin_center << ")";
    return os;
}