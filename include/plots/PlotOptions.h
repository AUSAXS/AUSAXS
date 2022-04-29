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

            void set(std::map<std::string, std::any> options);

        private: 
            enum class option {COLOR, ALPHA, MARKER_STYLE, LINE_WIDTH, MARKER_SIZE, DRAW_LINE, DRAW_ERROR, DRAW_MARKER, USE_EXISTING_AXES, TITLE, XLABEL, YLABEL};
            const inline static std::map<option, std::vector<std::string>> aliases = {
                {option::COLOR, {"color", "colour", "c"}},
                {option::ALPHA, {"alpha"}},
                {option::MARKER_STYLE, {"markerstyle", "marker_style", "ms"}},
                {option::LINE_WIDTH, {"linewidth", "line_width", "lw"}},
                {option::MARKER_SIZE, {"markersize", "marker_size", "s"}},
                {option::DRAW_LINE, {"drawline", "draw_line"}},
                {option::DRAW_ERROR, {"drawerror", "draw_error", "drawerrors", "draw_errors"}},
                {option::DRAW_MARKER, {"drawmarker", "draw_marker", "drawmarkers", "draw_markers"}},
                {option::USE_EXISTING_AXES, {"useexistingaxes", "use-existing-axes", "use_existing_axes"}},
                {option::TITLE, {"title"}},
                {option::XLABEL, {"xlabel"}},
                {option::YLABEL, {"ylabel"}}
            };

            void parse(std::string key, std::any val);

            void set(const option opt, std::string name, std::any val);
    };
    
    [[maybe_unused]]
    static void draw(const std::shared_ptr<TGraph> graph, const PlotOptions& options) {
        // prepare visuals
        graph->SetLineColorAlpha(options.color, options.alpha);
        graph->SetMarkerColorAlpha(options.color, options.alpha);
        graph->SetMarkerStyle(options.marker_style);
        graph->SetLineWidth(options.line_width);
        graph->SetMarkerSize(options.marker_size);

        if (!options.title.empty()) {graph->SetTitle(options.title.c_str());};
        if (!options.xlabel.empty()) {
            graph->GetXaxis()->SetTitle(options.xlabel.c_str());
            graph->GetXaxis()->CenterTitle();
        }
        if (!options.ylabel.empty()) {
            graph->GetYaxis()->SetTitle(options.ylabel.c_str());
            graph->GetYaxis()->CenterTitle();
            graph->GetYaxis()->SetTitleOffset(1.2);
        }

        // draw it
        std::string draw_options = options.use_existing_axes ? "SAME " : "A";
        if (options.draw_line) {draw_options += "L";}
        if (options.draw_markers) {draw_options += "P";}
        if (options.draw_errors) {graph->DrawClone(draw_options.c_str());}
        else {graph->TGraph::DrawClone(draw_options.c_str());}
    }

    [[maybe_unused]]
    static void draw(const std::shared_ptr<TH1D> hist, const PlotOptions& options) {
        // prepare visuals
        hist->SetLineColorAlpha(options.color, options.alpha);
        hist->SetMarkerColorAlpha(options.color, options.alpha);
        hist->SetMarkerStyle(options.marker_style);
        hist->SetLineWidth(options.line_width);
        hist->SetMarkerSize(options.marker_size);

        if (!options.title.empty()) {hist->SetTitle(options.title.c_str());};
        if (!options.xlabel.empty()) {
            hist->GetXaxis()->SetTitle(options.xlabel.c_str());
            hist->GetXaxis()->CenterTitle();
        }
        if (!options.ylabel.empty()) {
            hist->GetYaxis()->SetTitle(options.ylabel.c_str());
            hist->GetYaxis()->CenterTitle();
            hist->GetYaxis()->SetTitleOffset(1.2);
        }

        hist->DrawClone("HIST L");
    }
}