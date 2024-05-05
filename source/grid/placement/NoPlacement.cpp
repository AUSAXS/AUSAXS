/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/placement/NoPlacement.h>
#include <grid/GridMember.h>
#include <data/record/Water.h>

std::vector<grid::GridMember<data::record::Water>> grid::NoPlacement::place() const {
    return {};
}