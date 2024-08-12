#pragma once

#include <math/Vector.h>
#include <dataset/DatasetFwd.h>
#include <utility/Axis.h>
#include <utility/TypeTraits.h>

namespace hist {
    /**
     * @brief A representation of a histogram. 
     */
    class Histogram {
        public:
            Histogram() = default;
            Histogram(const Histogram& rhs) = default;
            Histogram(Histogram&& rhs) noexcept = default;
            Histogram& operator=(const Histogram& rhs) = default;
            Histogram& operator=(Histogram&& rhs) noexcept = default;
            virtual ~Histogram() = default;

            /**
             * @brief Construct a new histogram based on a list of bin values. 
             *        The axis will be initialized to [0, p.size()].
             * 
             * @param p The bin values. 
             */
            Histogram(const Vector<double>& p) noexcept;

            /**
             * @brief Construct a new histogram based on a list of bin values and the axis it spans. 
             * 
             * @param p The bin values. 
             * @param axis The axis they span. 
             */
            Histogram(const Vector<double>& p, const Axis& axis);

            /**
             * @brief Construct a new histogram based on a list of bin values and the axis it spans. 
             * 
             * @param p The bin values. 
             * @param axis The axis they span. 
             */
            Histogram(std::vector<double>&& p_tot, const Axis& axis);

            /**
             * @brief Construct a new empty histogram. 
             * 
             * @param axis The axis range. 
             */
            Histogram(const Axis& axis) noexcept;

            /**
             * @brief Reduce the view axis to show only the non-zero area. 
             *        Minimum size is 10 units.
             */
            void shorten_axis(unsigned int min_size = 10);

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
             * @brief Get the axis of this Histogram.
             */
            const Axis& get_axis() const;

            /**
             * @brief Get the total histogram.
             */
            virtual const std::vector<double>& get_counts() const;

            std::vector<double>& get_counts();

            /**
             * @brief Get the count at a specific bin.
             */
            const double& get_count(unsigned int i) const;

            /**
             * @brief Add a count to a specific bin.
             */
            void add_count(unsigned int i, double count);

            /**
             * @brief Set the count at a specific bin.
             */
            void set_count(unsigned int i, double count);

            /**
             * @brief Set the counts of this histogram.
             */
            void set_count(const std::vector<double>& counts);

            /**
             * @brief Bin a value.
             *        The bin containing the value will be incremented by 1.
             */
            void bin(double value);

            /**
             * @brief Bin a list of values.
             *        The bin containing each value will be incremented by 1.
             */
            void bin(const std::vector<double>& values);

            /**
             * @brief Get the count at a specific bin.
             */
            double& get_count(unsigned int i);

            const double& index(unsigned int i) const; // @copydoc get_count(unsigned int i) const

            double& index(unsigned int i); // @copydoc get_count(unsigned int i)

            /**
             * @brief Get the spanned range of this histogram. 
             */
            [[nodiscard]] Limit span_y() const noexcept;

            /**
             * @brief Get the positive spanned range of this histogram.
             *        This can be useful for setting log ranges. 
             */
            [[nodiscard]] Limit span_y_positive() const noexcept;

            /**
             * @brief Get the size of this Histogram.
             */
            [[nodiscard]] unsigned int size() const noexcept;

            [[nodiscard]] virtual std::string to_string() const noexcept;

            [[nodiscard]] SimpleDataset as_dataset() const;

            /**
             * @brief Normalize the histogram to have a sum of 1. 
             */
            void normalize();

            Histogram& operator+=(const Histogram& rhs);
            Histogram& operator-=(const Histogram& rhs);
            Histogram& operator*=(double rhs);
            double& operator[](int i);
            double operator[](int i) const;
            bool operator==(const Histogram& rhs) const;

        protected:
            mutable Vector<double> p;   // The bin values. 
            Axis axis;                  // The axis spanned by this histogram. 
    };

    Histogram operator*(const Histogram& lhs, double rhs);
    Histogram operator-(const Histogram& lhs, const Histogram& rhs);
    Histogram operator+(const Histogram& lhs, const Histogram& rhs);

    static_assert(supports_nothrow_move_v<Histogram>, "Histogram should be noexcept move constructible.");
}