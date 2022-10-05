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
            int marker_style = 8;           // Marker style
            int line_style = 1;             // Line style
            unsigned int line_width = 1;    // Line width
            double marker_size = 0.5;       // Marker size
            bool draw_line = true;          // Draw a line through the points
            bool draw_errors = false;       // Draw error bars if possible
            bool draw_markers = false;      // Draw a marker for each point
            bool draw_bars = false;         // Draw bars for a histogram.
            bool use_existing_axes = false; // Draw with existing axes. Must be false for the first plot on each canvas. 
            bool logx = false;              // Log scale for the x-axis. Only valid if use_existing_axes is false.  
            bool logy = false;              // Log scale for the y-axis. Only valid if use_existing_axes is false. 
            Limit ylimits;                  // Limits on the y-axis
            Limit xlimits;                  // Limits on the x-axis

            // cosmetic
            std::string title = "";         // Title
            std::string xlabel = "";        // Label for the x-axis
            std::string ylabel = "";        // Label for the y-axis
            double xlabel_offset = 0;       // Offset for the x-axis label
            double ylabel_offset = 0;       // Offset for the y-axis label

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
                make_shared({"useexistingaxes", "use-existing-axes", "use_existing_axes", "share_axes", "share_axis"}, use_existing_axes),
                make_shared({"title"}, title),
                make_shared({"xlabel"}, xlabel),
                make_shared({"ylabel"}, ylabel),
                make_shared({"logx", "log_x"}, logx),
                make_shared({"logy", "log_y"}, logy),
                make_shared({"xlim", "x_lim", "xlimits", "xlimit"}, xlimits),
                make_shared({"ylim", "y_lim", "ylimits", "ylimit"}, ylimits),
                make_shared({"xlabeloffset", "x_label_offset", "xlabel_offset", "x_offset"}, xlabel_offset),
                make_shared({"ylabeloffset", "y_label_offset", "ylabel_offset", "y_offset"}, ylabel_offset),
            };

            void parse(std::string key, std::any val);
    };    

    /**
     * @brief A small wrapper class for PlotOptions, intended for inheritance. 
     */
    class PlotOptionWrapper {
        public:         
            /**
             * @brief Set the plot options for this dataset. 
             */
            void set_plot_options(const plots::PlotOptions& options);

            /**
             * @brief Add plot options for this dataset.
             */
            void add_plot_options(const std::map<std::string, std::any>& options);

            /**
             * @brief Add plot options for this dataset, forcing the specified style. 
             *        Accepted styles: "line", "marker", "errors".
             */
            void add_plot_options(std::string style, std::map<std::string, std::any> options = {});

            /**
             * @brief Add plot options for this dataset, forcing the specified color. 
             */
            void add_plot_options(int color, std::map<std::string, std::any> options = {});

            /**
             * @brief Set the plot color for this dataset. 
             */
            void set_plot_color(int color);
        
            /**
             * @brief Get the current plot options.
             */
            plots::PlotOptions get_plot_options() const;

        protected: 
            plots::PlotOptions options;
    };
}