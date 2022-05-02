#pragma once

#include <string>
#include <any>

#include <TROOT.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TAxis.h>

namespace plots {
    /**
     * @brief Defines the plot style used for a single plot.
     */
    class PlotOptions {
        public:
            PlotOptions() {}

            PlotOptions(int color);

            /**
             * @brief Create a new set of plot settings. 
             * 
             * @param style Should be either "markers" or "line", depending on which style is required. 
             * @param options The remaining options. If both markers and lines are needed, set both to true here. 
             */
            PlotOptions(std::string style, std::map<std::string, std::any> options);

            int color = kBlack;             // Color
            double alpha = 1;               // Opacity
            int marker_style = 7;           // Marker style
            unsigned int line_width = 1;    // Line width
            unsigned int marker_size = 1;   // Marker size
            bool draw_line = false;         // Draw a line through the points
            bool draw_errors = true;        // Draw error bars if possible
            bool draw_markers = false;      // Draw a marker for each point
            bool use_existing_axes = false; // Draw with existing axes. Must be false for the first plot on each canvas. 
            std::string title = "";         // Title
            std::string xlabel = "";        // Label for the x-axis
            std::string ylabel = "";        // Label for the y-axis

            PlotOptions& set(std::map<std::string, std::any> options);

            PlotOptions& set(std::string style, std::map<std::string, std::any> options);

        private: 
            enum class option {COLOR, ALPHA, MARKER_STYLE, LINE_WIDTH, MARKER_SIZE, DRAW_LINE, DRAW_ERROR, DRAW_MARKER, USE_EXISTING_AXES, TITLE, XLABEL, YLABEL};
            const inline static std::map<option, std::vector<std::string>> aliases = {
                {option::COLOR, {"color", "colour", "c"}},
                {option::ALPHA, {"alpha"}},
                {option::MARKER_STYLE, {"markerstyle", "marker_style", "ms"}},
                {option::LINE_WIDTH, {"linewidth", "line_width", "lw"}},
                {option::MARKER_SIZE, {"markersize", "marker_size", "s"}},
                {option::DRAW_LINE, {"drawline", "draw_line", "drawlines", "draw_lines", "line", "lines"}},
                {option::DRAW_ERROR, {"drawerror", "draw_error", "drawerrors", "draw_errors", "error", "errors"}},
                {option::DRAW_MARKER, {"drawmarker", "draw_marker", "drawmarkers", "draw_markers", "draw_points", "drawpoints", "markers", "marker"}},
                {option::USE_EXISTING_AXES, {"useexistingaxes", "use-existing-axes", "use_existing_axes", "share_axes", "share_axis"}},
                {option::TITLE, {"title"}},
                {option::XLABEL, {"xlabel"}},
                {option::YLABEL, {"ylabel"}}
            };

            void parse(std::string key, std::any val);

            void set(const option opt, std::string name, std::any val);
    };    
}