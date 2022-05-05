#include <plots/Plot.h>
#include <utility/Exceptions.h>

#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>

void plots::draw(const std::shared_ptr<TGraph> graph, const PlotOptions& options, const std::shared_ptr<TCanvas> canvas) {
    // handle colors and markers
    graph->SetLineColorAlpha(options.color, options.alpha);
    graph->SetMarkerColorAlpha(options.color, options.alpha);
    graph->SetMarkerStyle(options.marker_style);
    graph->SetLineWidth(options.line_width);
    graph->SetMarkerSize(options.marker_size);

    // handle title & labels
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

    // handle log scale
    if (options.logx) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogx();
    }
    if (options.logy) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogy();
    }

    // prepare draw options
    std::string draw_options = options.use_existing_axes ? "SAME " : "A";
    if (options.draw_line) {draw_options += "L";}
    if (options.draw_markers) {draw_options += "P";}

    // draw it
    if (options.draw_errors) {graph->DrawClone(draw_options.c_str());}
    else {graph->TGraph::DrawClone(draw_options.c_str());}
}

void plots::draw(const Dataset& data, const std::shared_ptr<TCanvas> canvas) {
    draw(data, data.plot_options, canvas);
}

void plots::draw(const Dataset& data, const PlotOptions& options, const std::shared_ptr<TCanvas> canvas) {
    std::shared_ptr<TGraph> graph;
    if (!data.yerr.empty()) {graph = std::make_shared<TGraphErrors>(data.x.size(), data.x.data(), data.y.data(), data.xerr.data(), data.yerr.data());}
    else {graph = std::make_shared<TGraph>(data.x.size(), data.x.data(), data.y.data());}
    draw(graph, options, canvas);
}

void plots::draw(const std::shared_ptr<TH1D> hist, const PlotOptions& options, const std::shared_ptr<TCanvas> canvas) {
    // handle colors and markers
    hist->SetLineColorAlpha(options.color, options.alpha);
    hist->SetMarkerColorAlpha(options.color, options.alpha);
    hist->SetMarkerStyle(options.marker_style);
    hist->SetLineWidth(options.line_width);
    hist->SetMarkerSize(options.marker_size);

    // handle title & labels
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

    // handle log scale
    if (options.logx) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogx();
    }
    if (options.logy) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogy();
    }

    hist->DrawClone("HIST L");
}