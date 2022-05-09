#include <plots/PlotOptions.h>
#include <utility/Exceptions.h>

#include <algorithm>
#include <typeindex>

#include <TGraphErrors.h>

using namespace plots;
using std::string;

PlotOptions::PlotOptions() : draw_line(true) {}

PlotOptions::PlotOptions(int color) : color(color), draw_line(true) {}

PlotOptions::PlotOptions(const PlotOptions& opt) {*this = opt;}

PlotOptions::PlotOptions(string style, std::map<std::string, std::any> options) {
    parse(style, true);
    set(options);
}

PlotOptions& PlotOptions::set(string style, std::map<string, std::any> options) {
    draw_line = draw_markers = draw_errors = false;
    parse(style, true);
    set(options);
    return *this;
}

PlotOptions& PlotOptions::set(std::map<string, std::any> options) {
    std::for_each(options.begin(), options.end(), [this] (const auto& opt) {parse(opt.first, opt.second);});
    return *this;
}

PlotOptions& PlotOptions::set(int color, std::map<std::string, std::any> options) {
    parse("color", color);
    set(options);
    return *this;
}

void PlotOptions::parse(string key, std::any val) {
    for (const auto& opt : options) {
        for (const auto& alias : opt->aliases) {
            if (alias == key) {
                opt->parse(val);
                return;
            }
        }
    }
    throw except::parse_error("Error in PlotOptions::parse: Unknown option \"" + key + "\".");
}

template<>
void PlotOptions::SmartOption<string>::parse(const std::any val) {
    if (std::type_index{typeid(const char*)} == val.type()) {
        value = string(std::any_cast<const char*>(val));
    } else if (std::type_index{typeid(string)} == val.type()) {
        value = std::any_cast<string>(val);
    } else {
        throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be a string. Received \"" + std::string(typeid(val.type()).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<bool>::parse(const std::any val) {
    if (std::type_index{typeid(bool)} == val.type()) {
        value = std::any_cast<bool>(val);
    } else {
        throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be a boolean. Received \"" + std::string(typeid(val.type()).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<int>::parse(const std::any val) {
    if (std::type_index{typeid(int)} == val.type()) {
        value = std::any_cast<int>(val);
    } else if (std::type_index{typeid(EColor)} == val.type()) {
        value = std::any_cast<EColor>(val);
    } else {
        throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be an integer. Received \"" + std::string(typeid(val.type()).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<unsigned int>::parse(const std::any val) {
    if (std::type_index{typeid(unsigned int)} == val.type()) {
        value = std::any_cast<unsigned int>(val);
    } else {
        throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be an integer. Received \"" + std::string(typeid(val.type()).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<double>::parse(const std::any val) {
    if (std::type_index{typeid(double)} == val.type()) {
        value = std::any_cast<double>(val);
    } else {
        throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be a double. Received \"" + std::string(typeid(val.type()).name()) + "\".");
    }
}

PlotOptions& PlotOptions::operator=(const PlotOptions& opt) {
    color = opt.color;
    alpha = opt.alpha; 
    marker_style = opt.marker_style; 
    line_width = opt.line_width; 
    marker_size = opt.marker_size; 
    draw_line = opt.draw_line;
    draw_errors = opt.draw_errors;
    draw_markers = opt.draw_markers;
    use_existing_axes = opt.use_existing_axes;
    logx = opt.logx; 
    logy = opt.logy;
    title = opt.title; 
    xlabel = opt.xlabel; 
    ylabel = opt.ylabel; 
    return *this;
}