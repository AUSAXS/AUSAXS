#pragma once

#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>

namespace ausaxs::grid {
    template<typename C>
    concept grid_member_t = std::is_base_of_v<data::AtomFF, C> || std::is_base_of_v<data::Water, C>;

	template<grid_member_t T> class GridMember;
	class PlacementStrategy;
	class CullingStrategy;
}