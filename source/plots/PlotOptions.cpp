#include <plots/PlotOptions.h>
#include <utility/Exceptions.h>

#include <algorithm>
#include <typeindex>
#include <initializer_list>

using namespace plots;
using std::string;

double inf = std::numeric_limits<double>::infinity();

PlotOptions::PlotOptions() : draw_line(true) {}

PlotOptions::PlotOptions(const PlotOptions& opt) {*this = opt;}

PlotOptions::PlotOptions(string style, std::map<std::string, std::any> options) {
    draw_line = draw_markers = draw_errors = false;
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
void PlotOptions::SmartOption<Limit>::parse(const std::any val) {
    // handle actual Limit case: Limit(1.5, 2.5)
    if (std::type_index{typeid(Limit)} == val.type()) {
        value = std::any_cast<Limit>(val);
    } 

    // handle double list initializer like {1.5, 2.5}
    else if (std::type_index{typeid(std::vector<double>)} == val.type()) {
        std::vector<double> vals(std::any_cast<std::vector<double>>(val));
        if (vals.size() != 2) {throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must contain two values. Received \"" + std::to_string(vals.size()) + "\".");}
        value = Limit(vals[0], vals[1]);
    } 
    
    // handle integer list initializer like {1, 2}
    else if (std::type_index{typeid(std::vector<int>)} == val.type()) {
        std::vector<int> vals(std::any_cast<std::vector<int>>(val));
        if (vals.size() != 2) {throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must contain two values. Received \"" + std::to_string(vals.size()) + "\".");}
        value = Limit(vals[0], vals[1]);
    } 
    
    // otherwise throw
    else { throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be a pair of two values. Received \"" + std::string(typeid(val.type()).name()) + "\".");}
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
    } else {
        throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be an integer. Received \"" + std::string(typeid(val.type()).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<unsigned int>::parse(const std::any val) {
    if (std::type_index{typeid(unsigned int)} == val.type()) {
        value = std::any_cast<unsigned int>(val);
    } else if (std::type_index{typeid(int)} == val.type()) {
        int parsed_val = std::any_cast<int>(val);
        if (parsed_val < 0) {throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be strictly positive.");}
        value = parsed_val;
    } else {
        throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be an integer. Received \"" + std::string(typeid(val.type()).name()) + "\".");
    }
}

template<>
void PlotOptions::SmartOption<double>::parse(const std::any val) {
    if (std::type_index{typeid(double)} == val.type()) {
        value = std::any_cast<double>(val);
    } else if (std::type_index{typeid(int)} == val.type()) {
        value = std::any_cast<int>(val);
    } else {
        throw except::invalid_argument("Error in PlotOptions::set: Option \"" + aliases[0] + "\" must be a double. Received \"" + std::string(typeid(val.type()).name()) + "\".");
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
    use_existing_axes = opt.use_existing_axes;
    logx = opt.logx; 
    logy = opt.logy;
    title = opt.title; 
    xlabel = opt.xlabel; 
    ylabel = opt.ylabel; 
    draw_bars = opt.draw_bars;
    ylimits = opt.ylimits;
    xlimits = opt.xlimits;
    xlabel_offset = opt.xlabel_offset;
    ylabel_offset = opt.ylabel_offset;

    return *this;
}

std::string PlotOptions::to_string() const {
    std::stringstream ss;
    ss << "PlotOptions " << std::endl;
    ss << "    color "          << color << std::endl;
    ss << "    marker_style "   << marker_style << std::endl;
    ss << "    line_style "     << line_style << std::endl;
    ss << "    line_width "     << line_width << std::endl;
    ss << "    marker_size "    << marker_size << std::endl;
    ss << "    draw_line "      << draw_line << std::endl;
    ss << "    draw_errors "    << draw_errors << std::endl;
    ss << "    draw_markers "   << draw_markers << std::endl;
    ss << "    logx "           << logx << std::endl;
    ss << "    logy "           << logy << std::endl;
    ss << "    title "          << title << std::endl;
    ss << "    xlabel "         << xlabel << std::endl;
    ss << "    ylabel "         << ylabel << std::endl;
    ss << "    ylimits "        << ylimits.min << ylimits.max << std::endl;
    ss << "    xlimits "        << xlimits.min << ylimits.max << std::endl;
    return ss.str();
}

void PlotOptionWrapper::set_plot_options(const plots::PlotOptions& options) {this->options = options;}

void PlotOptionWrapper::add_plot_options(std::map<std::string, std::any>& options) {this->options.set(options);}

void PlotOptionWrapper::add_plot_options(std::string style, std::map<std::string, std::any> options) {this->options.set(style, options);}

void PlotOptionWrapper::set_plot_color(std::string color) {this->options.color = color;}

plots::PlotOptions PlotOptionWrapper::get_plot_options() const {return options;}