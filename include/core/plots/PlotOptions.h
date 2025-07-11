// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/Limit.h>
#include <plots/Styles.h>

#include <vector>
#include <string>
#include <any>
#include <unordered_map>
#include <memory>

namespace ausaxs::plots {
    extern double inf;

    struct option {
        // style
        static inline std::string color = "color";
        static inline std::string alpha = "alpha";
        static inline std::string line_style = "linestyle";
        static inline std::string marker_style = "markerstyle";
        static inline std::string line_width = "linewidth";
        static inline std::string marker_size = "markersize";

        // draw options
        static inline std::string draw_line = "line";
        static inline std::string draw_errors = "error";
        static inline std::string draw_markers = "marker";
        static inline std::string draw_bars = "bars";
        static inline std::string draw_fill = "fill";
        static inline std::string zorder = "zorder";

        // labels
        static inline std::string title = "title";
        static inline std::string xlabel = "xlabel";
        static inline std::string ylabel = "ylabel";
        static inline std::string zlabel = "zlabel";
        static inline std::string x2label = "x2label";
        static inline std::string legend = "legend";

        // axis settings
        static inline std::string logx = "logx";
        static inline std::string logy = "logy";
        static inline std::string xlimits = "xlim";
        static inline std::string ylimits = "ylim";
        static inline std::string normalize = "normalize";
    };

    /**
     * @brief Defines the plot style used for a single plot.
     */
    class PlotOptions {
        public:
            PlotOptions();
            PlotOptions(const PlotOptions& opt);

            /**
             * @brief Create a new set of plot settings. 
             * 
             * @param options The remaining options. If both markers and lines are needed, set both to true here. 
             */
            PlotOptions(std::unordered_map<std::string, std::any> options);

            /**
             * @brief Create a new set of plot settings. 
             * 
             * @param style Should be either "markers" or "line", depending on which style is required. 
             * @param options The remaining options. If both markers and lines are needed, set both to true here. 
             */
            PlotOptions(const style::DrawStyle& style, std::unordered_map<std::string, std::any> options);

            PlotOptions(std::initializer_list<std::pair<std::string, std::any>> options);

            PlotOptions& set(std::unordered_map<std::string, std::any> options);

            PlotOptions& set(const style::DrawStyle& style, std::unordered_map<std::string, std::any> options = {});

            PlotOptions& operator=(const PlotOptions& opt);

            std::string to_string() const;

            // remember to add new options to the equality operator overload
            style::Color color = "k";               // Color. Default is black.
            double alpha = 1;                       // Opacity
            std::string marker_style = ".";         // Marker style
            style::LineStyle line_style = "-";      // Line style
            double line_width = 1;                  // Line width
            double marker_size = 5;                 // Marker size
            bool draw_line = true;                  // Draw a line through the points
            bool draw_errors = false;               // Draw error bars if possible
            bool draw_markers = false;              // Draw a marker for each point
            bool draw_bars = false;                 // Draw bars for a histogram.
            bool logx = false;                      // Log scale for the x-axis. Only valid if use_existing_axes is false.  
            bool logy = false;                      // Log scale for the y-axis. Only valid if use_existing_axes is false. 
            bool normalize = false;                 // Normalize the y-axis to 1.
            Limit ylimits;                          // Limits on the y-axis
            Limit xlimits;                          // Limits on the x-axis
            int zorder = 0;                         // Z-order for the plot. Higher zorder plots are drawn on top of lower zorder plots.

            // cosmetic
            std::string title = "";                 // Title
            std::string xlabel = "x";               // Label for the x-axis
            std::string x2label = "";               // Label for the secondary x-axis
            std::string ylabel = "y";               // Label for the y-axis
            std::string zlabel = "z";               // Label for the z-axis
            std::string legend = "";                // Legend entry

        private: 
            struct ISmartOption {
                ISmartOption(const std::vector<std::string>& aliases) : aliases(aliases) {}
                virtual ~ISmartOption() = default;

                virtual void parse(const std::any& val) = 0;

                std::vector<std::string> aliases;
            };

            template<typename T>
            struct SmartOption : ISmartOption {
                SmartOption(const std::vector<std::string>& aliases, T& value) : ISmartOption(aliases), value(value) {}
                ~SmartOption() override = default;

                void parse(const std::any& val) override;

                T& value;
            };

            template<typename T>
            std::shared_ptr<SmartOption<T>> make_shared(std::vector<std::string> aliases, T& val) {
                return std::make_shared<SmartOption<T>>(aliases, val);
            }

            const std::vector<std::shared_ptr<ISmartOption>> options = {
                make_shared({option::color, "colour", "c"}, color),
                make_shared({option::alpha}, alpha),
                make_shared({option::line_style, "line_style", "ls"}, line_style),
                make_shared({option::marker_style, "marker_style", "ms"}, marker_style),
                make_shared({option::line_width, "line_width", "lw"}, line_width),
                make_shared({option::marker_size, "marker_size", "s"}, marker_size),
                make_shared({option::draw_line, "lines"}, draw_line),
                make_shared({option::draw_errors, "errors"}, draw_errors),
                make_shared({option::draw_markers, "markers", "point", "points"}, draw_markers),
                make_shared({option::draw_bars}, draw_bars),
                make_shared({option::title}, title),
                make_shared({option::xlabel}, xlabel),
                make_shared({option::ylabel}, ylabel),
                make_shared({option::zlabel}, zlabel),
                make_shared({option::logx, "log_x"}, logx),
                make_shared({option::logy, "log_y"}, logy),
                make_shared({option::normalize, "normalize"}, normalize),
                make_shared({option::xlimits, "x_lim", "xlimits", "xlimit", "xlims", "xrange"}, xlimits),
                make_shared({option::ylimits, "y_lim", "ylimits", "ylimit", "ylims", "yrange"}, ylimits),
                make_shared({option::legend}, legend),
                make_shared({option::zorder}, zorder)
            };

            void parse(const std::string& key, std::any val);
    };
}