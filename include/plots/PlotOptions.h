#pragma once

#include <vector>
#include <string>
#include <any>
#include <map>
#include <memory>

#include <utility/Axis.h>

extern double inf;

namespace plots {
    /**
     * @brief Defines the plot style used for a single plot.
     */
    class PlotOptions {
        public:
            PlotOptions();

            PlotOptions(int color);

            /**
             * @brief Create a new set of plot settings. 
             * 
             * @param style Should be either "markers" or "line", depending on which style is required. 
             * @param options The remaining options. If both markers and lines are needed, set both to true here. 
             */
            PlotOptions(std::string style, std::map<std::string, std::any> options);

            /**
             * @brief Copy constructor.
             */
            PlotOptions(const PlotOptions& opt);

            PlotOptions& set(std::map<std::string, std::any> options);

            PlotOptions& set(std::string style, std::map<std::string, std::any> options = {});

            PlotOptions& set(int color, std::map<std::string, std::any> options = {});

            PlotOptions& operator=(const PlotOptions& opt);

            // remember to add new options to the equality operator overload
            int color = 1;                  // Color. Default is kBlack = 1
            double alpha = 1;               // Opacity
            int marker_style = 7;           // Marker style
            unsigned int line_width = 1;    // Line width
            double marker_size = 1;         // Marker size
            bool draw_line = true;          // Draw a line through the points
            bool draw_errors = false;       // Draw error bars if possible
            bool draw_markers = false;      // Draw a marker for each point
            bool draw_bars = false;         // Draw bars for a histogram.
            bool use_existing_axes = false; // Draw with existing axes. Must be false for the first plot on each canvas. 
            bool logx = false;              // Log scale for the x-axis. Only valid if use_existing_axes is false.  
            bool logy = false;              // Log scale for the y-axis. Only valid if use_existing_axes is false. 
            Limit ylimits;                  // Limits on the y-axis
            Limit xlimits;                  // Limits on the x-axis
            std::string title = "";         // Title
            std::string xlabel = "";        // Label for the x-axis
            std::string ylabel = "";        // Label for the y-axis

        private: 
            enum class option {COLOR, ALPHA, MARKER_STYLE, LINE_WIDTH, MARKER_SIZE, DRAW_LINE, DRAW_ERROR, DRAW_MARKER, DRAW_BARS, USE_EXISTING_AXES, TITLE, XLABEL, YLABEL, LOGX, LOGY, YLIM, XLIM};

            struct ISmartOption {
                ISmartOption(option opt, std::vector<std::string> aliases) : opt(opt), aliases(aliases) {}
                virtual ~ISmartOption() = default;

                virtual void parse(const std::any val) = 0;

                option opt;
                std::vector<std::string> aliases;
            };

            template<typename T>
            struct SmartOption : ISmartOption {
                SmartOption(option opt, std::vector<std::string> aliases, T& value) : ISmartOption(opt, aliases), value(value) {}
                ~SmartOption() override = default;

                void parse(const std::any val) override;

                T& value;
            };

            template<typename T>
            std::shared_ptr<SmartOption<T>> make_shared(option opt, std::vector<std::string> aliases, T& val) {
                return std::make_shared<SmartOption<T>>(opt, aliases, val);
            }

            const std::vector<std::shared_ptr<ISmartOption>> options = {
                make_shared(option::COLOR, {"color", "colour", "c"}, color),
                make_shared(option::ALPHA, {"alpha"}, alpha),
                make_shared(option::MARKER_STYLE, {"markerstyle", "marker_style", "ms"}, marker_style),
                make_shared(option::LINE_WIDTH, {"linewidth", "line_width", "lw"}, line_width),
                make_shared(option::MARKER_SIZE, {"markersize", "marker_size", "s"}, marker_size),
                make_shared(option::DRAW_LINE, {"line", "lines"}, draw_line),
                make_shared(option::DRAW_ERROR, {"error", "errors"}, draw_errors),
                make_shared(option::DRAW_MARKER, {"marker", "markers", "point", "points"}, draw_markers),
                make_shared(option::DRAW_BARS, {"bars", "bars"}, draw_bars),
                make_shared(option::USE_EXISTING_AXES, {"useexistingaxes", "use-existing-axes", "use_existing_axes", "share_axes", "share_axis"}, use_existing_axes),
                make_shared(option::TITLE, {"title"}, title),
                make_shared(option::XLABEL, {"xlabel"}, xlabel),
                make_shared(option::YLABEL, {"ylabel"}, ylabel),
                make_shared(option::LOGX, {"logx", "log_x"}, logx),
                make_shared(option::LOGY, {"logy", "log_y"}, logy),
                make_shared(option::XLIM, {"xlim", "x_lim"}, xlimits),
                make_shared(option::YLIM, {"ylim", "y_lim"}, ylimits)
            };

            void parse(std::string key, std::any val);
    };    
}