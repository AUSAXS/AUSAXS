#include <plots/PlotOptions.h>
#include <utility/Exceptions.h>

#include <algorithm>
#include <typeindex>

using namespace plots;
using std::string;

PlotOptions::PlotOptions(int color) : color(color) {}

PlotOptions& PlotOptions::set(std::map<string, std::any> options) {
    std::for_each(options.begin(), options.end(), [this] (const auto& opt) {parse(opt.first, opt.second);});
    return *this;
}

void PlotOptions::parse(string key, std::any val) {
    for (const auto& e : aliases) {
        for (const auto& alias : e.second) {
            if (alias == key) {
                set(e.first, alias, val);
                return;
            }
        }
    }
    throw except::parse_error("Error in PlotOptions::parse: Unknown option \"" + key + "\".");
}

void PlotOptions::set(const option opt, string name, std::any val) {
    static auto get_string = [] (string name, std::any val) {
        if (std::type_index{typeid(const char*)} == val.type()) {
            return string(std::any_cast<const char*>(val));
        } else if (std::type_index{typeid(string)} == val.type()) {
            return std::any_cast<string>(val);
        } else {
            throw except::invalid_argument("Error in PlotOptions::set: Option \"" + name + "\" must be a string. Received \"" + std::string(typeid(val.type()).name()) + "\".");
        }
    };

    static auto get_int = [] (string name, std::any val) {
        if (std::type_index{typeid(int)} == val.type()) {
            return std::any_cast<int>(val);
        } else {
            throw except::invalid_argument("Error in PlotOptions::set: Option \"" + name + "\" must be an integer. Received \"" + std::string(typeid(val.type()).name()) + "\".");
        }
    };

    static auto get_double = [] (string name, std::any val) {
        if (std::type_index{typeid(double)} == val.type()) {
            return std::any_cast<double>(val);
        } else {
            throw except::invalid_argument("Error in PlotOptions::set: Option \"" + name + "\" must be a double. Received \"" + std::string(typeid(val.type()).name()) + "\".");
        }
    };

    static auto get_bool = [] (string name, std::any val) {
        if (std::type_index{typeid(bool)} == val.type()) {
            return std::any_cast<bool>(val);
        } else {
            throw except::invalid_argument("Error in PlotOptions::set: Option \"" + name + "\" must be a boolean. Received \"" + std::string(typeid(val.type()).name()) + "\".");
        }
    };

    switch (opt) {
        case option::COLOR: // these can both be just integers and ROOT EColors
            if (std::type_index{typeid(int)} == val.type()) {
                color = std::any_cast<int>(val);
            } else if (std::type_index{typeid(EColor)} == val.type()) {
                color = std::any_cast<EColor>(val);
            } else {
                throw except::invalid_argument("Error in PlotOptions::set: Color option must be an integer. Received \"" + std::string(typeid(val.type()).name()) + "\".");
            }
            break;
        case option::ALPHA: 
            alpha = get_double(name, val);
            break;
        case option::MARKER_STYLE: 
            marker_style = get_int(name, val);
            break;
        case option::LINE_WIDTH: 
            line_width = get_int(name, val);
            break;
        case option::MARKER_SIZE: 
            marker_size = get_int(name, val);
            break;
        case option::DRAW_LINE: 
            draw_line = get_bool(name, val);
            break;
        case option::DRAW_ERROR: 
            draw_errors = get_bool(name, val);
            break;
        case option::DRAW_MARKER: 
            draw_markers = get_bool(name, val);
            break;
        case option::USE_EXISTING_AXES: 
            use_existing_axes = get_bool(name, val);
            break;
        case option::TITLE: 
            title = get_string(name, val);
            break;
        case option::XLABEL: 
            xlabel = get_string(name, val);
            break;
        case option::YLABEL: 
            ylabel = get_string(name, val);
            break;
        default: 
            throw except::unexpected("Error in PlotOptions::set: Received unhandled valid option \"" + name + "\".");
    }
}

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