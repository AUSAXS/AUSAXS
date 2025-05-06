#pragma once

#include <grid/expansion/GridExpander.h>

namespace ausaxs::grid::expander {
    struct FullExpander : public GridExpander {
        FullExpander(observer_ptr<Grid> grid);

        void expand_volume(GridMember<data::AtomFF>& atom);
        void expand_volume(GridMember<data::Water>& atom);

        void deflate_volume(GridMember<data::AtomFF>& atom);
        void deflate_volume(GridMember<data::Water>& atom);
    };
}