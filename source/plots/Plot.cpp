#include <plots/Plot.h>
#include <utility/Exceptions.h>
#include <histogram/Histogram.h>

#include <TPad.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>

void plots::detail::handle_log(const PlotOptions& options, const std::shared_ptr<TCanvas> canvas) {
    if (options.logx && !options.use_existing_axes) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogx();
    }
    if (options.logy && !options.use_existing_axes) {
        if (canvas == nullptr) {throw except::nullptr_error("Error in Plot::draw: Can only set log scale if canvas is provided.");}
        canvas->SetLogy();
    }
}

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

    // handle xlimits
    if (!options.xlimits.empty()) {
        Limit limits = options.xlimits;
        // replace infs
        if (std::isinf(limits.min)) {limits.min = graph->GetXaxis()->GetXmin();}
        if (std::isinf(limits.max)) {limits.max = graph->GetXaxis()->GetXmax();}

        // set new range
        graph->GetXaxis()->SetRangeUser(options.xlimits.min, options.xlimits.max);
    }

    // handle ylimits
    if (!options.ylimits.empty()) {
        if (!std::isinf(options.ylimits.min)) {graph->SetMinimum(options.ylimits.min);} 
        if (!std::isinf(options.ylimits.max)) {graph->SetMaximum(options.ylimits.max);}
    }

    // handle log scale
    detail::handle_log(options, canvas);

    // prepare draw options
    std::string draw_options = options.use_existing_axes ? "SAME " : "A";
    if (options.draw_line) {draw_options += "L";}
    if (options.draw_markers || options.draw_errors) {draw_options += "P";}
    if (options.draw_bars) {throw except::invalid_argument("Error in Plot::Draw: Invalid option for Dataset: \"bars\"");}

    // draw it
    graph->DrawClone(draw_options.c_str());
}

std::shared_ptr<TGraph> plots::detail::graph(const Dataset& data, const plots::PlotOptions& options) {
    std::shared_ptr<TGraph> graph;
    if (options.draw_errors) {
        graph = std::make_shared<TGraphErrors>(data.size(), data.x().to_vector().data(), data.y().to_vector().data(), data.xerr().to_vector().data(), data.yerr().to_vector().data());
    } else {
        graph = std::make_shared<TGraph>(data.size(), data.x().to_vector().data(), data.y().to_vector().data());
    } 
    return graph;    
}

std::shared_ptr<TGraph> plots::detail::graph(const Dataset& data) {
    return graph(data, data.get_plot_options());
}

std::shared_ptr<TGraph> plots::detail::graph(const SimpleDataset& data, const plots::PlotOptions& options) {
    std::shared_ptr<TGraph> graph;
    if (options.draw_errors) {
        graph = std::make_shared<TGraphErrors>(data.size(), data.x().to_vector().data(), data.y().to_vector().data(), nullptr, data.yerr().to_vector().data());
    } else {
        graph = std::make_shared<TGraph>(data.size(), data.x().to_vector().data(), data.y().to_vector().data());
    }
    return graph;
}

std::shared_ptr<TGraph> plots::detail::graph(const SimpleDataset& data) {
    return graph(data, data.get_plot_options());
}

void plots::draw(const Multiset& data, const std::shared_ptr<TCanvas> canvas) {
    PlotOptions options = data[0].get_plot_options();

    TMultiGraph graphs;
    for (const Dataset& d : data) {
        PlotOptions options = d.get_plot_options();

        std::string draw_options;
        if (options.draw_line) {draw_options += "L";}
        if (options.draw_markers || options.draw_errors) {draw_options += "P";}
        if (options.draw_bars) {throw except::invalid_argument("Error in Plot::Draw: Invalid option for Dataset: \"bars\"");}
        TGraph* g = new TGraph(*detail::graph(d, options)); // TMultiGraph owns its TGraph pointers

        g->SetLineColorAlpha(options.color, options.alpha);
        g->SetMarkerColorAlpha(options.color, options.alpha);
        g->SetMarkerStyle(options.marker_style);
        g->SetLineWidth(options.line_width);
        g->SetMarkerSize(options.marker_size);

        graphs.Add(g, draw_options.c_str());
    }    

    std::string labels = options.title + ";" + options.xlabel + ";" + options.ylabel;

    // set title & labels
    graphs.SetTitle(labels.c_str());
    graphs.GetXaxis()->CenterTitle();
    graphs.GetYaxis()->CenterTitle();

    // handle xlimits
    if (!options.xlimits.empty()) {
        Limit limits = options.xlimits;
        // replace infs
        if (std::isinf(limits.min)) {limits.min = graphs.GetXaxis()->GetXmin();}
        if (std::isinf(limits.max)) {limits.max = graphs.GetXaxis()->GetXmax();}

        // set new range
        graphs.GetXaxis()->SetRangeUser(options.xlimits.min, options.xlimits.max);
    }

    // handle ylimits
    if (!options.ylimits.empty()) {
        if (!std::isinf(options.ylimits.min)) {graphs.SetMinimum(options.ylimits.min);} 
        if (!std::isinf(options.ylimits.max)) {graphs.SetMaximum(options.ylimits.max);}
    }

    // handle log
    detail::handle_log(options, canvas);

    graphs.DrawClone("A");
}

void plots::draw(const Dataset& data, const std::shared_ptr<TCanvas> canvas) {
    draw(data, data.get_plot_options(), canvas);
}

void plots::draw(const SimpleDataset& data, const std::shared_ptr<TCanvas> canvas) {
    draw(data, data.get_plot_options(), canvas);
}

void plots::draw(const SimpleDataset& data, const PlotOptions& options, const std::shared_ptr<TCanvas> canvas) {
    auto g = detail::graph(data, options);
    draw(g, options, canvas);
}

void plots::draw(const Dataset& data, const PlotOptions& options, const std::shared_ptr<TCanvas> canvas) {
    auto g = detail::graph(data, options);
    draw(g, options, canvas);
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

    // handle xlimits
    if (!options.xlimits.empty()) {
        Limit limits = options.xlimits;
        // replace infs
        if (std::isinf(limits.min)) {limits.min = hist->GetXaxis()->GetXmin();}
        if (std::isinf(limits.max)) {limits.max = hist->GetXaxis()->GetXmax();}

        // set new range
        hist->GetXaxis()->SetRangeUser(options.xlimits.min, options.xlimits.max);
    }

    // handle ylimits
    if (!options.ylimits.empty()) {
        if (!std::isinf(options.ylimits.min)) {hist->SetMinimum(options.ylimits.min);} 
        if (!std::isinf(options.ylimits.max)) {hist->SetMaximum(options.ylimits.max);}
    }

    // handle log scale
    detail::handle_log(options, canvas);

    // prepare draw options
    std::string draw_options = options.use_existing_axes ? "HIST SAME " : "HIST ";
    if (options.draw_bars) {draw_options += "B";}
    if (options.draw_line) {draw_options += "L";}
    if (options.draw_markers || options.draw_errors) {draw_options += "P";}
    if (options.draw_bars) {throw except::invalid_argument("Error in Plot::Draw: Invalid option for Dataset: \"bars\"");}

    hist->DrawClone("HIST L");
}