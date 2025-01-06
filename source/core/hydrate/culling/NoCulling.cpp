/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/NoCulling.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>

void ausaxs::hydrate::NoCulling::cull(std::span<grid::GridMember<data::Water>>&) const {return;}