#pragma once

#include <utility/Axis.h>
#include <math/Vector.h>
#include <plots/PlotOptions.h>

namespace hist {
    /**
     * @brief \class Histogram. 
     * 
     * A representation of a histogram. 
     */
    class Histogram {
        public:
            /**
             * @brief Default constructor.
             */
            Histogram() noexcept {}

            /**
             * @brief Constructor. 
             * 
             * Construct a new histogram based on a list of bin values. 
             * Note that the axis will not be initialized. 
             * 
             * @param p The bin values. 
             */
            Histogram(const Vector<double>& p) noexcept;

            /**
             * @brief Constructor.
             * 
             * Construct a new histogram based on a list of bin values and the axis it spans. 
             * 
             * @param p The bin values. 
             * @param axis The axis they span. 
             */
            Histogram(const Vector<double>& p, const Axis& axis) noexcept;

            /**
             * @brief Constructor.
             * 
             * Construct a new empty histogram. 
             * 
             * @param axis The axis range. 
             */
            Histogram(const Axis& axis) noexcept;

            /**
             * @brief Add another Histogram to this one.
             */
            Histogram& operator+=(const Histogram& rhs) noexcept;

            /**
             * @brief Subtract another Histogram from this one.
             */
            Histogram& operator-=(const Histogram& rhs) noexcept;

            /**
             * @brief Reduce the view axis to show only the non-zero area. 
             *        Minimum size is 10 units.
             */
            void shorten_axis();

            /**
             * @brief Automatically generate an axis containing all elements. 
             */
            void generate_axis(unsigned int size = 100);

            /**
             * @brief Set the axis of this Histogram.
             */
            void set_axis(const Axis& axis) noexcept;

            /**
             * @brief Get the minimum value. 
             */
            [[nodiscard]] Limit span() const noexcept;

            /**
             * @brief Get the positive span of this histogram.
             *        This can be useful for setting log ranges. 
             */
            [[nodiscard]] Limit span_positive() const noexcept;

            /**
             * @brief Get the size of this Histogram.
             */
            [[nodiscard]] size_t size() const noexcept;

            /**
             * @brief Overwrite the plot options for this Histogram.
             */
            void set_plot_options(const plots::PlotOptions& options);

            /**
             * @brief Change the plot options for this Histogram.
             */
            void add_plot_options(const std::map<std::string, std::any>& options);

            /**
             * @brief Change the plot options for this Histogram.
             * 
             * @param style The plotting style. Should be one of the accepted variations of "markers" or "line". 
             * @param options The other plot options.
             */
            void add_plot_options(std::string style, std::map<std::string, std::any> options = {});

            /**
             * @brief Change the plot options for this dataset.
             * 
             * @param color The color.
             * @param options The other plot options.
             */
            void add_plot_options(int color, std::map<std::string, std::any> options = {});

            Vector<double> p;                // The bin values. 
            Axis axis;                       // The axis spanned by this histogram. 
            plots::PlotOptions plot_options; // The plot options used to draw this Histogram.
    };
}