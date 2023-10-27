#pragma once

#include <initializer_list>
#include <vector>
#include <array>
#include <string>

class Limit;

/**
 * @brief A representation of a one-dimensional axis. 
 * 		  The left and right limits are inclusive.
 */
class Axis {
	public:
		/**
		 * @brief Default constructor. 
		 */
		Axis() noexcept;

		/**
		 * @brief Construct a new axis based on a Limit and the number of bins. 
		 * 
		 * @param limits The limits on the axis. 
		 * @param bins The number of equidistant bins. 
		 */
		Axis(const Limit& limits, int bins) noexcept;

		/**
		 * @brief Construct a new axis based on explicit minimum and maximum values, along with the number of bins. 
		 * 
		 * @param xmin The minimum value spanned by this axis. 
		 * @param xmax The maximum value spanned by this axis. 
		 * @param bins The number of equidistant bins. 
		 */
		constexpr Axis(double xmin, double xmax, int bins) noexcept : bins(bins), min(xmin), max(xmax)  {}

		/**
		 * @brief List initializer. 
		 * 
		 * Initialize this Axis with an initializer_list. 
		 */
		Axis& operator=(std::initializer_list<double> list) noexcept;

		/**
		 * @brief Get a string representation of this object. 
		 */
		std::string to_string() const noexcept;

		/**
		 * @brief Equality operator. 
		 * 
		 * Check if this object is equal to another. 
		 */
		bool operator==(const Axis& rhs) const noexcept;

		/**
		 * @brief Inequality operator. 
		 * 
		 * Check if this object is different from another. 
		 */
		bool operator!=(const Axis& rhs) const noexcept;

		/**
		 * @brief Get the bin width. 
		 */
		constexpr double width() const noexcept {return (max-min)/bins;}

		/**
		 * @brief Get the span of this Axis.
		 */
		constexpr double span() const noexcept {return max-min;}

		/**
		 * @brief Get the bin width.
		 */
		constexpr double step() const noexcept {return width();}

		/**
		 * @brief Get the bin index for a given value.
		 * 		  Returns bins in the range [0, bins]
		 */
		unsigned int get_bin(double value) const noexcept;

		/**
		 * @brief Get the axis value for a given bin.
		 * 		  Returns values in the range [min, max]
		 */
		double get_bin_value(unsigned int bin) const noexcept;

		/**
		 * @brief Get a sub-axis of this Axis.
		 * 		  The closest match to the specified minimum and maximum values are used.
		 * 		  This guarantees that the returned Axis is an exact subset of this Axis.
		 */
		Axis sub_axis(double min, double max) const noexcept;

		/**
		 * @brief Resize this Axis to a new number of bins.
		 * 		  The maximum value is adjusted to keep the bin width constant.
		 */
		void resize(unsigned int bins) noexcept;

		/**
		 * @brief Get a vector representation of this Axis.
		 * 
		 * @param shift Specify the amount to shift each bin by. Using 0.5 will return the center values of each bin.
		 */
		std::vector<double> as_vector(double shift = 0) const noexcept {
			std::vector<double> v(bins);
			double w = width();
			double new_min = min + shift*w;
			for (unsigned int i = 0; i < bins; ++i) {
				v[i] = new_min + i*w;
			}
			return v;
		}

		/**
		 * @brief Get an array representation of this Axis. 
		 * 
		 * @param shift Specify the amount to shift each bin by. Using 0.5 will return the center values of each bin.
		 * @tparam size The size of the array.
		 */
		template<unsigned int size>
		constexpr std::array<double, size> as_array(double shift = 0) const noexcept {
			std::array<double, size> v;
			double w = width();
			double new_min = min + shift*w;
			for (unsigned int i = 0; i < bins; ++i) {
				v[i] = new_min + i*w;
			}
			return v;
		}

		/**
		 * @brief Check if this Axis is empty.
		 */
		[[nodiscard]] bool empty() const noexcept;

		/**
		 * @brief Get the limits of this Axis.
		 */
		Limit limits() const noexcept;

		unsigned int bins; 	// The number of equidistant bins. 
		double min;        	// The minimum value spanned by this Axis. 
		double max;        	// The maximum value spanned by this Axis. 
};