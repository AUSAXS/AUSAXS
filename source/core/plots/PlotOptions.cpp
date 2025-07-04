// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <plots/PlotOptions.h>
#include <utility/Exceptions.h>

#include <typeindex>
#include <algorithm>
#include <sstream>

using namespace ausaxs::plots;

double inf = std::numeric_limits<double>::infinity();

PlotOptions::PlotOptions() : draw_line(true) {}
PlotOptions::PlotOptions(const PlotOptions& opt) {*this = opt;}

PlotOptions::PlotOptions(std::unordered_map<std::string, std::any> options) : draw_line(true) {set(std::move(options));}

PlotOptions::PlotOptions(std::initializer_list<std::pair<std::string, std::any>> options) : draw_line(true) {
    std::unordered_map<std::string, std::any> opts;
    for (const auto& opt : options) {
        opts[opt.first] = opt.second;
    }
    set(std::move(opts));
}

PlotOptions::PlotOptions(const style::DrawStyle& style, std::unordered_map<std::string, std::any> options) {
    draw_line = draw_markers = draw_errors = false;
    parse(style, true);
    set(std::move(options));
}

PlotOptions& PlotOptions::set(const style::DrawStyle& style, std::unordered_map<std::string, std::any> options) {
    draw_line = draw_markers = draw_errors = false;
    parse(style, true);
    set(std::move(options));
    return *this;
}

PlotOptions& PlotOptions::set(std::unordered_map<std::string, std::any> options) {
    std::for_each(options.begin(), options.end(), [this] (const auto& opt) {parse(opt.first, opt.second);});
    return *this;
}

void PlotOptions::parse(const std::string& key, std::any val) {
    for (const auto& opt : options) {
        for (const auto& alias : opt->aliases) {
            if (alias == key) {
                opt->parse(val);
                return;
            }
        }
    }
    throw except::parse_error("PlotOptions::parse: Unknown option \"" + key + "\".");
}

template<>
void PlotOptions::SmartOption<ausaxs::Limit>::parse(const std::any& val) {
    // handle actual Limit case: Limit(1.5, 2.5)
    if (std::type_index{typeid(Limit)} == val.type()) {
        value = std::any_cast<Limit>(val);
    } 

    // handle double list initializer like {1.5, 2.5}
    else if (std::type_index{typeid(std::vector<double>)} == val.type()) {
        std::vector<double> vals(std::any_cast<std::vector<double>>(val));
        if (vals.size() != 2) {throw except::invalid_argument("PlotOptions::set: Option \"" + aliases[0] + "\" must contain two values. Received \"" + std::to_string(vals.size()) + "\".");}
        value = Limit(vals[0], vals[1]);
    } 
    
    // handle integer list initializer like {1, 2}
    else if (std::type_index{typeid(std::vector<int>)} == val.type()) {
        std::vector<int> vals(std::any_cast<std::vector<int>>(val));
        if (vals.size() != 2) {throw except::invalid_argument("PlotOptions::set: Option \"" + aliases[0] + "\" must contain two values. Received \"" + std::to_string(vals.size()) + "\".");}
        value = Limit(vals[0], vals[1]);
    } 
    
    // otherwise throw
    else {
        auto& t = val.type();
        throw except::invalid_argument("PlotOptions::set: Option \"" + aliases[0] + "\" must be a pair of two values. Received \"" + std::string(typeid(t).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<std::string>::parse(const std::any& val) {
    if (std::type_index{typeid(const char*)} == val.type()) {
        value = std::string(std::any_cast<const char*>(val));
    } else if (std::type_index{typeid(std::string)} == val.type()) {
        value = std::any_cast<std::string>(val);
    } else if (std::type_index{typeid(style::Color)} == val.type()) {
        value = std::any_cast<style::Color>(val);
    } else if (std::type_index{typeid(style::LineStyle)} == val.type()) {
        value = std::any_cast<style::LineStyle>(val);
    } else if (std::type_index{typeid(style::DrawStyle)} == val.type()) {
        value = std::any_cast<style::DrawStyle>(val);
    } else {
        auto& t = val.type();
        throw except::invalid_argument("PlotOptions::set: Option \"" + aliases[0] + "\" must be a string. Received \"" + std::string(typeid(t).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<bool>::parse(const std::any& val) {
    if (std::type_index{typeid(bool)} == val.type()) {
        value = std::any_cast<bool>(val);
    } else {
        auto& t = val.type();
        throw except::invalid_argument("PlotOptions::set: Option \"" + aliases[0] + "\" must be a boolean. Received \"" + std::string(typeid(t).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<int>::parse(const std::any& val) {
    if (std::type_index{typeid(int)} == val.type()) {
        value = std::any_cast<int>(val);
    } else {
        auto& t = val.type();
        throw except::invalid_argument("PlotOptions::set: Option \"" + aliases[0] + "\" must be an integer. Received \"" + std::string(typeid(t).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<unsigned int>::parse(const std::any& val) {
    if (std::type_index{typeid(unsigned int)} == val.type()) {
        value = std::any_cast<unsigned int>(val);
    } else if (std::type_index{typeid(int)} == val.type()) {
        int parsed_val = std::any_cast<int>(val);
        if (parsed_val < 0) {throw except::invalid_argument("PlotOptions::set: Option \"" + aliases[0] + "\" must be strictly positive.");}
        value = parsed_val;
    } else {
        auto& t = val.type();
        throw except::invalid_argument("PlotOptions::set: Option \"" + aliases[0] + "\" must be an integer. Received \"" + std::string(typeid(t).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<double>::parse(const std::any& val) {
    if (std::type_index{typeid(double)} == val.type()) {
        value = std::any_cast<double>(val);
    } else if (std::type_index{typeid(int)} == val.type()) {
        value = std::any_cast<int>(val);
    } else {
        auto& t = val.type();
        throw except::invalid_argument("PlotOptions::set: Option \"" + aliases[0] + "\" must be a double. Received \"" + std::string(typeid(t).name()) + "\".");
    }
}

PlotOptions& PlotOptions::operator=(const PlotOptions& opt) {
    color = opt.color;
    alpha = opt.alpha; 
    marker_style = opt.marker_style; 
    line_style = opt.line_style;
    line_width = opt.line_width; 
    marker_size = opt.marker_size; 
    draw_line = opt.draw_line;
    draw_errors = opt.draw_errors;
    draw_markers = opt.draw_markers;
    logx = opt.logx; 
    logy = opt.logy;
    title = opt.title; 
    xlabel = opt.xlabel; 
    x2label = opt.x2label;
    ylabel = opt.ylabel; 
    zlabel = opt.zlabel;
    draw_bars = opt.draw_bars;
    ylimits = opt.ylimits;
    xlimits = opt.xlimits;
    legend = opt.legend;
    zorder = opt.zorder;
    normalize = opt.normalize;
    return *this;
}

std::string PlotOptions::to_string() const {
    std::stringstream ss;
    ss << "PlotOptions " << "\n"
        << "color "          << color << "\n"
        << "marker_style "   << marker_style << "\n"
        << "line_style "     << line_style << "\n"
        << "line_width "     << line_width << "\n"
        << "marker_size "    << marker_size << "\n"
        << "draw_line "      << draw_line << "\n"
        << "draw_errors "    << draw_errors << "\n"
        << "draw_markers "   << draw_markers << "\n"
        << "logx "           << logx << "\n"
        << "logy "           << logy << "\n"
        << "title "          << title << "\n"
        << "xlabel "         << xlabel << "\n"
        << "x2label "        << x2label << "\n"
        << "ylabel "         << ylabel << "\n"
        << "zlabel "         << zlabel << "\n"
        << "legend "         << legend << "\n"
        << "normalize "      << normalize << "\n"
        << "ylimits "        << ylimits.min << " " << ylimits.max << "\n"
        << "xlimits "        << xlimits.min << " " << xlimits.max << "\n"
        << "zorder "         << zorder << std::endl;
    return ss.str();
}