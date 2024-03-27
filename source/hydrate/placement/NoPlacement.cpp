/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/placement/NoPlacement.h>
#include <hydrate/GridMember.h>
#include <data/record/Water.h>

std::vector<grid::GridMember<data::record::Water>> grid::NoPlacement::place() const {
    return {};
}