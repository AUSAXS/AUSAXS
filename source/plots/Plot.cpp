#include <plots/Plot.h>

#include <TGraphErrors.h>

void plots::draw(const std::shared_ptr<TGraph> graph, const PlotOptions& options) {
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

void plots::draw(const Dataset& data) {
    draw(data, data.plot_options);
}

void plots::draw(const Dataset& data, const PlotOptions& options) {
    std::shared_ptr<TGraph> graph;
    if (!data.yerr.empty()) {graph = std::make_shared<TGraphErrors>(data.x.size(), data.x.data(), data.y.data(), data.xerr.data(), data.yerr.data());}
    else {graph = std::make_shared<TGraph>(data.x.size(), data.x.data(), data.y.data());}
    draw(graph, options);
}

void plots::draw(const std::shared_ptr<TH1D> hist, const PlotOptions& options) {
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