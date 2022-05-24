#include <plots/Plot.h>
#include <utility/Exceptions.h>
#include <histogram/Histogram.h>

#include <TPad.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>

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
    if (options.logx && !options.use_existing_axes) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogx();
    }
    if (options.logy && !options.use_existing_axes) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogy();
    }

    // handle axis limits
    if (!options.xlimits.empty()) {
        Limit limits = options.xlimits;
        // check if new bounds are within current bounds
        if (!(graph->GetXaxis()->GetXmin() < limits.min && limits.max < graph->GetXaxis()->GetXmax())) {
            // replace infs
            if (std::isinf(options.xlimits.min)) {limits.min = graph->GetXaxis()->GetXmin();}
            if (std::isinf(options.xlimits.max)) {limits.max = graph->GetXaxis()->GetXmax();}

            // set new range
            graph->GetXaxis()->SetRangeUser(limits.min, limits.max);
        }
    }
    if (!options.ylimits.empty()) {
        Limit limits = options.ylimits;
        // check if new bounds are within current bounds
        if (!(graph->GetYaxis()->GetXmin() < limits.min && limits.max < graph->GetYaxis()->GetXmax())) {
            // set range
            if (!std::isinf(options.ylimits.min)) {graph->GetHistogram()->SetMinimum(limits.min);} 
            if (!std::isinf(options.ylimits.max)) {graph->GetHistogram()->SetMaximum(limits.max);}
        }
    }

    // prepare draw options
    std::string draw_options = options.use_existing_axes ? "SAME " : "A";
    if (options.draw_line) {draw_options += "L";}
    if (options.draw_markers) {draw_options += "P";}
    if (options.draw_bars) {throw except::invalid_argument("Error in Plot::Draw: Invalid option for Dataset: \"bars\"");}

    // draw it
    if (options.draw_errors) {graph->DrawClone(draw_options.c_str());}
    else {graph->TGraph::DrawClone(draw_options.c_str());}
}

void plots::draw(const Multiset& data, const std::shared_ptr<TCanvas> canvas) {
    TMultiGraph graph;
    for (const Dataset& d : data) {
        PlotOptions options = d.plot_options;

        std::string draw_options;
        if (options.draw_line) {draw_options += "L";}
        if (options.draw_markers) {draw_options += "P";}
        if (options.draw_bars) {throw except::invalid_argument("Error in Plot::Draw: Invalid option for Dataset: \"bars\"");}
        TGraph* g = new TGraph(*d.plot()); // TMultiGraph owns its TGraph pointers

        graph.Add(g, draw_options.c_str());
    }    

    graph.DrawClone("A");
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

void plots::draw(const hist::Histogram& hist, const PlotOptions& options, const std::shared_ptr<TCanvas> canvas) {
    if (hist.axis.empty()) {throw except::invalid_argument("Error in Plot::draw: The axis of a histogram must be defined before plotting.");}
    std::shared_ptr<TH1D> h = std::make_unique<TH1D>("hist", "hist", hist.axis.bins, hist.axis.min, hist.axis.max);
    std::for_each(hist.p.begin(), hist.p.end(), [&h] (double val) {h->Fill(val);});
    draw(h, options, canvas);
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

    // handle axis limits
    if (!options.xlimits.empty()) {
        Limit limits;
        if (std::isinf(options.xlimits.min)) {limits.min = hist->GetXaxis()->GetXmin();}
        if (std::isinf(options.xlimits.max)) {limits.max = hist->GetXaxis()->GetXmax();}
        hist->GetXaxis()->SetRangeUser(limits.min, limits.max);
    }
    if (!options.ylimits.empty()) {
        Limit limits;
        if (std::isinf(options.ylimits.min)) {limits.min = hist->GetYaxis()->GetXmin();}
        if (std::isinf(options.ylimits.max)) {limits.max = hist->GetYaxis()->GetXmax();}
        hist->GetXaxis()->SetRangeUser(limits.min, limits.max);
    }

    // handle log scale
    if (options.logx && !options.use_existing_axes) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogx();
    }
    if (options.logy && !options.use_existing_axes) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogy();
    }

    // prepare draw options
    std::string draw_options = options.use_existing_axes ? "HIST SAME " : "HIST ";
    if (options.draw_bars) {draw_options += "B";}
    if (options.draw_line) {draw_options += "L";}
    if (options.draw_markers) {draw_options += "P";}
    if (options.draw_bars) {throw except::invalid_argument("Error in Plot::Draw: Invalid option for Dataset: \"bars\"");}

    hist->DrawClone("HIST L");
}