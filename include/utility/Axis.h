#pragma once

#include <utility/Limit.h>

#include <initializer_list>
#include <vector>
#include <string>

/**
 * @brief A representation of an axis. 
 */
class Axis {
	public:
		/**
		 * @brief Default constructor. 
		 */
		Axis() noexcept;

		/**
		 * @brief Constructor. 
		 * 
		 * Construct a new axis based on a Limit and the number of bins. 
		 * 
		 * @param limits The limits on the axis. 
		 * @param bins The number of equidistant bins. 
		 */
		Axis(const Limit& limits, int bins) noexcept;

		/**
		 * @brief Constructor. 
		 * 
		 * Construct a new axis based on explicit minimum and maximum values, along with the number of bins. 
		 * 
		 * @param bins The number of equidistant bins. 
		 * @param xmin The minimum value spanned by this axis. 
		 * @param xmax The maximum value spanned by this axis. 
		 */
		Axis(int bins, double xmin, double xmax) noexcept;

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
		double width() const noexcept;

		/**
		 * @brief Get the span of this Axis.
		 */
		double span() const noexcept;

		/**
		 * @brief Get the bin width.
		 */
		double step() const noexcept;

		void resize(unsigned int bins) noexcept;

		/**
		 * @brief Get a vector representation of this Axis.
		 */
		std::vector<double> as_vector() const noexcept;

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