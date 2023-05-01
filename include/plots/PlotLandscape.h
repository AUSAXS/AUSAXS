#pragma once

#include <plots/Plot.h>
#include <mini/detail/Landscape.h>

namespace plots {
    /**
     * @brief Plot a landscape.
     */
    class PlotLandscape : public Plot {
        public:
            /**
             * @brief Constructor.
             * 
             * @param data The landscape which will be plotted. 
             * @param path The path to the folder where the plot will be saved. 
             */
            PlotLandscape(const mini::Landscape& data, const io::File& path);

            /**
             * @brief Destructor. 
             */
            ~PlotLandscape() override;
        
            /**
             * @brief Plot and save the input dataset at the specified location. 
             * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
             */
            static void quick_plot(const mini::Landscape& data, const io::File& path);
    };
}