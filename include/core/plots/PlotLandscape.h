// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <plots/Plot.h>
#include <mini/MiniFwd.h>

namespace ausaxs::plots {
    /**
     * @brief Plot a specific \class Landscape object.
     */
    class PlotLandscape : public Plot {
        public:
            /**
             * @brief Constructor.
             * 
             * @param data The landscape which will be plotted. 
             * @param path The path to the folder where the plot will be saved. 
             */
            PlotLandscape(const mini::Landscape& data, const PlotOptions& options, const io::File& path);

            /**
             * @brief Destructor. 
             */
            ~PlotLandscape() override;
        
            /**
             * @brief Plot and save the input dataset at the specified location. 
             * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
             */
            static void quick_plot(const mini::Landscape& data, const PlotOptions& options, const io::File& path);
    };
}