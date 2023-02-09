#pragma once

#include <vector>
#include <string>
#include <any>
#include <map>
#include <memory>

#include <utility/Axis.h>

extern double inf;

struct style {
    typedef std::string LineStyle;
    typedef std::string DrawStyle;
	typedef std::string Color;

    struct color {
        inline static Color black = "k";
        inline static Color blue = "tab:blue";
        inline static Color orange = "tab:orange";
        inline static Color green = "tab:green";
        inline static Color red = "tab:red";
        inline static Color purple = "tab:purple";
        inline static Color brown = "tab:brown";
        inline static Color pink = "tab:pink";
        inline static Color gray = "tab:gray";
        inline static Color olive = "tab:olive";
        inline static Color cyan = "tab:cyan";
    };

    struct line {
        inline static LineStyle solid = "-";
        inline static LineStyle dashed = "--";
        inline static LineStyle dotted = ":";
        inline static LineStyle dashdot = "-.";
    };

    struct draw {
        inline static DrawStyle line = "line";
        inline static DrawStyle hist = "hist";
        inline static DrawStyle points = "points";
        inline static DrawStyle errors = "errors";
    };
};

namespace plots {
    /**
     * @brief Defines the plot style used for a single plot.
     */
    class PlotOptions {
        public:
            PlotOptions();

            /**
             * @brief Create a new set of plot settings. 
             * 
             * @param style Should be either "markers" or "line", depending on which style is required. 
             * @param options The remaining options. If both markers and lines are needed, set both to true here. 
             */
            PlotOptions(style::DrawStyle style, std::map<std::string, std::any> options);

            /**
             * @brief Create a new set of plot settings. 
             * 
             * @param style Should be either "markers" or "line", depending on which style is required. 
             * @param options The remaining options. If both markers and lines are needed, set both to true here. 
             */
            PlotOptions(style, std::map<std::string, std::any> options);

            /**
             * @brief Copy constructor.
             */
            PlotOptions(const PlotOptions& opt);

            PlotOptions& set(std::map<std::string, std::any> options);

            PlotOptions& set(style::DrawStyle style, std::map<std::string, std::any> options = {});

            PlotOptions& operator=(const PlotOptions& opt);

            std::string to_string() const;

            // remember to add new options to the equality operator overload
            style::Color color = "k";               // Color. Default is black.
            double alpha = 1;                       // Opacity
            std::string marker_style = ".";         // Marker style
            style::LineStyle line_style = "-";      // Line style
            unsigned int line_width = 1;            // Line width
            double marker_size = 5;                 // Marker size
            bool draw_line = true;                  // Draw a line through the points
            bool draw_errors = false;               // Draw error bars if possible
            bool draw_markers = false;              // Draw a marker for each point
            bool draw_bars = false;                 // Draw bars for a histogram.
            bool logx = false;                      // Log scale for the x-axis. Only valid if use_existing_axes is false.  
            bool logy = false;                      // Log scale for the y-axis. Only valid if use_existing_axes is false. 
            Limit ylimits;                          // Limits on the y-axis
            Limit xlimits;                          // Limits on the x-axis

            // cosmetic
            std::string title = "";         // Title
            std::string xlabel = "x";        // Label for the x-axis
            std::string ylabel = "y";        // Label for the y-axis
            std::string zlabel = "z";        // Label for the z-axis
            std::string legend = "";        // Legend entry

        private: 
            struct ISmartOption {
                ISmartOption(std::vector<std::string> aliases) : aliases(aliases) {}
                virtual ~ISmartOption() = default;

                virtual void parse(const std::any val) = 0;

                std::vector<std::string> aliases;
            };

            template<typename T>
            struct SmartOption : ISmartOption {
                SmartOption(std::vector<std::string> aliases, T& value) : ISmartOption(aliases), value(value) {}
                ~SmartOption() override = default;

                void parse(const std::any val) override;

                T& value;
            };

            template<typename T>
            std::shared_ptr<SmartOption<T>> make_shared(std::vector<std::string> aliases, T& val) {
                return std::make_shared<SmartOption<T>>(aliases, val);
            }

            const std::vector<std::shared_ptr<ISmartOption>> options = {
                make_shared({"color", "colour", "c"}, color),
                make_shared({"alpha"}, alpha),
                make_shared({"linestyle", "line_style", "ls"}, line_style),
                make_shared({"markerstyle", "marker_style", "ms"}, marker_style),
                make_shared({"linewidth", "line_width", "lw"}, line_width),
                make_shared({"markersize", "marker_size", "s"}, marker_size),
                make_shared({"line", "lines"}, draw_line),
                make_shared({"error", "errors"}, draw_errors),
                make_shared({"marker", "markers", "point", "points"}, draw_markers),
                make_shared({"bars", "bars"}, draw_bars),
                make_shared({"title"}, title),
                make_shared({"xlabel"}, xlabel),
                make_shared({"ylabel"}, ylabel),
                make_shared({"zlabel"}, zlabel),
                make_shared({"logx", "log_x"}, logx),
                make_shared({"logy", "log_y"}, logy),
                make_shared({"xlim", "x_lim", "xlimits", "xlimit"}, xlimits),
                make_shared({"ylim", "y_lim", "ylimits", "ylimit"}, ylimits),
                make_shared({"legend"}, legend)
            };

            void parse(std::string key, std::any val);
    };    

    /**
     * @brief A small wrapper class for PlotOptions, intended for inheritance. 
     */
    class Plottable {
        public:         
            /**
             * @brief Set the plot options for this dataset. 
             */
            void set_plot_options(const plots::PlotOptions& options);

            /**
             * @brief Add plot options for this dataset.
             */
            void add_plot_options(std::map<std::string, std::any> options);

            /**
             * @brief Add plot options for this dataset, forcing the specified style. 
             *        Accepted styles: "line", "marker", "errors".
             */
            void add_plot_options(std::string style, std::map<std::string, std::any> options = {});

            /**
             * @brief Set the plot color for this dataset. 
             */
            void set_plot_color(std::string color);
        
            /**
             * @brief Get the current plot options.
             */
            plots::PlotOptions get_plot_options() const;

        protected: 
            plots::PlotOptions options;
    };
}