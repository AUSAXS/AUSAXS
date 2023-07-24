#pragma once

#include <math/Vector.h>
#include <utility/Axis.h>
#include <plots/PlotOptions.h>

class SimpleDataset;
namespace hist {
    /**
     * @brief A representation of a histogram. 
     */
    class Histogram : public plots::Plottable {
        public:
            /**
             * @brief Default constructor.
             */
            Histogram() noexcept = default;

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

            virtual ~Histogram();

            /**
             * @brief Reduce the view axis to show only the non-zero area. 
             *        Minimum size is 10 units.
             */
            void shorten_axis(unsigned int min_size = 10);

            /**
             * @brief Extend the view axis to the given maximum value. 
             */
            void extend_axis(double qmax);

            /**
             * @brief Resize the number of bins in this histogram, keeping the width constant.
             */
            void resize(unsigned int bins);

            /**
             * @brief Automatically generate an axis containing all elements. 
             */
            void generate_axis();

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
            [[nodiscard]] unsigned int size() const noexcept;

            [[nodiscard]] virtual std::string to_string() const noexcept;

            [[nodiscard]] SimpleDataset as_dataset() const;

            Histogram& operator+=(const Histogram& rhs);
            Histogram& operator-=(const Histogram& rhs);
            Histogram& operator*=(double rhs);
            double& operator[](int i);
            double operator[](int i) const;
            bool operator==(const Histogram& rhs) const;

            Vector<double> p;   // The bin values. 
            Axis axis;          // The axis spanned by this histogram. 
    };

    Histogram operator*(const Histogram& lhs, double rhs);
    Histogram operator-(const Histogram& lhs, const Histogram& rhs);
    Histogram operator+(const Histogram& lhs, const Histogram& rhs);
}