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
    class Histogram : public plots::PlotOptionWrapper {
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

            Histogram& operator+=(const Histogram& rhs);
            Histogram& operator-=(const Histogram& rhs);
            double& operator[](const int i);
            double operator[](const int i) const;

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
             * @brief Get the spanned range of this histogram. 
             */
            [[nodiscard]] Limit span() const noexcept;

            /**
             * @brief Get the positive spanned range of this histogram.
             *        This can be useful for setting log ranges. 
             */
            [[nodiscard]] Limit span_positive() const noexcept;

            /**
             * @brief Get the size of this Histogram.
             */
            [[nodiscard]] size_t size() const noexcept;

            Vector<double> p;                // The bin values. 
            Axis axis;                       // The axis spanned by this histogram. 
    };
}