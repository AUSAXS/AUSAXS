#pragma once

#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <grid/detail/GridMember.h>
#include <grid/GridFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::grid::expander {
    struct GridExpander {
        GridExpander(observer_ptr<grid::Grid> grid);
        virtual ~GridExpander();

        /** 
         * @brief Expand a single member atom into an actual sphere.
         * 		  Only expands atoms if they have not already been expanded. 
         * 		  Complexity: O(1).
         */
        virtual void expand_volume(GridMember<data::AtomFF>& atom) = 0;

        /** 
         * @brief Expand a single member water molecule into an actual sphere.
         * 		  Only expands molecules if they have not already been expanded.
         * 		  Complexity: O(1).
         */
        virtual void expand_volume(GridMember<data::Water>& atom) = 0;

        /** 
         * @brief Deflate a single member atom into an actual sphere.
         * 		  Only deflates atoms if they have been expanded.
         * 		  Complexity: O(1).
         */
        virtual void deflate_volume(GridMember<data::AtomFF>& atom) = 0;

        /** 
         * @brief Deflate a single member atom into an actual sphere.
         * 		  Only deflates atoms if they have been expanded.
         * 		  Complexity: O(1).
         */
        virtual void deflate_volume(GridMember<data::Water>& atom) = 0;

        double to_x(int i) const;
        double to_y(int j) const;
        double to_z(int k) const;
        int& get_volume() const;

        observer_ptr<grid::Grid> grid;
    };
}