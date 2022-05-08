#pragma once

#include <data/Axis.h>
#include <math/Vector.h>

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
            Histogram() {}

            /**
             * @brief Constructor. 
             * 
             * Construct a new histogram based on a list of bin values. 
             * Note that the axis will not be initialized. 
             * 
             * @param p The bin values. 
             */
            Histogram(const Vector<double>& p);

            /**
             * @brief Constructor.
             * 
             * Construct a new histogram based on a list of bin values and the axis it spans. 
             * 
             * @param p The bin values. 
             * @param axis The axis they span. 
             */
            Histogram(const Vector<double>& p, const Axis& axis);

            /**
             * @brief Constructor.
             * 
             * Construct a new empty histogram. 
             * 
             * @param axis The axis range. 
             */
            Histogram(const Axis& axis);

            /**
             * @brief Add another Histogram to this one.
             */
            Histogram& operator+=(const Histogram& rhs);

            /**
             * @brief Subtract another Histogram from this one.
             */
            Histogram& operator-=(const Histogram& rhs);

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
            void set_axis(const Axis& axis);

            /**
             * @brief Get the size of this Histogram.
             */
            size_t size() const;

            Vector<double> p; // The bin values. 
            Axis axis;        // The axis spanned by this histogram. 
    };
}