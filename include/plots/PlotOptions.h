#pragma once

#include <string>
#include <initializer_list>

#include <TROOT.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TAxis.h>

#include <Dataset.h>

namespace plots {
    /**
     * @brief Defines the plot style used for a single plot.
     */
    class PlotOptions {
        public:
            PlotOptions() {}

            PlotOptions(int color) : color(color) {}

            int color = kBlack;
            double alpha = 1;
            int marker_style = 7;
            unsigned int linewidth = 1;
            unsigned int markersize = 1;
            bool draw_line = true;
            bool draw_markers = false;
            bool use_existing_axes = false;
            std::string title = "";
            std::string xlabel = "";
            std::string ylabel = "";
    };
    
    [[maybe_unused]]
    static void draw(const std::shared_ptr<TGraph> graph, const PlotOptions& options) {
        // prepare visuals
        graph->SetLineColorAlpha(options.color, options.alpha);
        graph->SetMarkerColorAlpha(options.color, options.alpha);
        graph->SetMarkerStyle(options.marker_style);
        graph->SetLineWidth(options.linewidth);
        graph->SetMarkerSize(options.markersize);

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
        if (options.draw_line) {graph->DrawClone(std::string(draw_options + "L").c_str());}
        if (options.draw_markers) {graph->DrawClone(std::string(draw_options + "P").c_str());}
    }

    [[maybe_unused]]
    static void draw(const std::shared_ptr<TH1D> hist, const PlotOptions& options) {
        // prepare visuals
        hist->SetLineColorAlpha(options.color, options.alpha);
        hist->SetMarkerColorAlpha(options.color, options.alpha);
        hist->SetMarkerStyle(options.marker_style);
        hist->SetLineWidth(options.linewidth);
        hist->SetMarkerSize(options.markersize);

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