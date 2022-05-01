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

            int color = kBlack;             // Color
            double alpha = 1;               // Opacity
            int marker_style = 7;           // Marker style
            unsigned int line_width = 1;     // Line width
            unsigned int marker_size = 1;    // Marker size
            bool draw_line = true;          // Draw a line through the points
            bool draw_errors = true;        // Draw error bars if possible
            bool draw_markers = false;      // Draw a marker for each point
            bool use_existing_axes = false; // Draw with existing axes. Must be false for the first plot on each canvas. 
            std::string title = "";         // Title
            std::string xlabel = "";        // Label for the x-axis
            std::string ylabel = "";        // Label for the y-axis

            PlotOptions& set(std::map<std::string, std::any> options);

        private: 
            enum class option {COLOR, ALPHA, MARKER_STYLE, LINE_WIDTH, MARKER_SIZE, DRAW_LINE, DRAW_ERROR, DRAW_MARKER, USE_EXISTING_AXES, TITLE, XLABEL, YLABEL};
            const inline static std::map<option, std::vector<std::string>> aliases = {
                {option::COLOR, {"color", "colour", "c"}},
                {option::ALPHA, {"alpha"}},
                {option::MARKER_STYLE, {"markerstyle", "marker_style", "ms"}},
                {option::LINE_WIDTH, {"linewidth", "line_width", "lw"}},
                {option::MARKER_SIZE, {"markersize", "marker_size", "s"}},
                {option::DRAW_LINE, {"drawline", "draw_line", "drawlines", "draw_lines"}},
                {option::DRAW_ERROR, {"drawerror", "draw_error", "drawerrors", "draw_errors"}},
                {option::DRAW_MARKER, {"drawmarker", "draw_marker", "drawmarkers", "draw_markers", "draw_points", "drawpoints"}},
                {option::USE_EXISTING_AXES, {"useexistingaxes", "use-existing-axes", "use_existing_axes"}},
                {option::TITLE, {"title"}},
                {option::XLABEL, {"xlabel"}},
                {option::YLABEL, {"ylabel"}}
            };

            void parse(std::string key, std::any val);

            void set(const option opt, std::string name, std::any val);
    };
    
    [[maybe_unused]]
    void draw(const std::shared_ptr<TGraph> graph, const PlotOptions& options);

    [[maybe_unused]]
    void draw(const std::shared_ptr<TGraph> graph);

    [[maybe_unused]]
    void draw(const std::shared_ptr<TH1D> hist, const PlotOptions& options);
}