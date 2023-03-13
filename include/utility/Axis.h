#pragma once

#include <utility/Limit.h>
#include <math/Vector3.h>

#include <vector>
#include <array>
#include <initializer_list>
#include <string>
#include <ostream>

/**
 * @brief A representation of an axis. 
 */
class Axis {
	public:
		/**
		 * @brief Default constructor. 
		 */
		Axis() noexcept : bins(0), min(0), max(0) {}

		/**
		 * @brief Constructor. 
		 * 
		 * Construct a new axis based on a Limit and the number of bins. 
		 * 
		 * @param limits The limits on the axis. 
		 * @param bins The number of equidistant bins. 
		 */
		Axis(const Limit& limits, int bins) noexcept : bins(bins), min(limits.min), max(limits.max) {}

		/**
		 * @brief Constructor. 
		 * 
		 * Construct a new axis based on explicit minimum and maximum values, along with the number of bins. 
		 * 
		 * @param bins The number of equidistant bins. 
		 * @param xmin The minimum value spanned by this axis. 
		 * @param xmax The maximum value spanned by this axis. 
		 */
		Axis(int bins, double xmin, double xmax) noexcept : bins(bins), min(xmin), max(xmax)  {}

		/**
		 * @brief List initializer. 
		 * 
		 * Initialize this Axis with an initializer_list. 
		 */
		Axis& operator=(std::initializer_list<double> list) noexcept {
			std::vector<double> d = list;
			bins = std::round(d[0]); 
			min = d[1];
			max = d[2];
			return *this;
		}

		/**
		 * @brief Get a string representation of this object. 
		 */
		std::string to_string() const noexcept {
			return "Axis: (" + std::to_string(min) + ", " + std::to_string(max) + ") with " + std::to_string(bins) + " bins";
		}

		/**
		 * @brief Equality operator. 
		 * 
		 * Check if this object is equal to another. 
		 */
		bool operator==(const Axis& rhs) const noexcept {
			if (bins != rhs.bins) {return false;}
			if (min != rhs.min) {return false;}
			if (max != rhs.max) {return false;}
			return true;
		}

		/**
		 * @brief Inequality operator. 
		 * 
		 * Check if this object is different from another. 
		 */
		bool operator!=(const Axis& rhs) const noexcept {return !operator==(rhs);}

		/**
		 * @brief Stream output operator. 
		 * 
		 * Allows this object to easily be output to a given stream. 
		 */
		friend std::ostream& operator<<(std::ostream& os, const Axis& axis) noexcept {os << axis.to_string(); return os;}

		/**
		 * @brief Get the bin width. 
		 */
		double width() const noexcept {return (max-min)/bins;}

		/**
		 * @brief Get the span of this Axis.
		 */
		double span() const noexcept {return max-min;}

		/**
		 * @brief Get the bin width.
		 */
		double step() const noexcept {return width();}

		/**
		 * @brief Get a vector representation of this Axis.
		 */
		std::vector<double> as_vector() const noexcept {
			std::vector<double> v(bins);
			double w = width();
			for (unsigned int i = 0; i < bins; i++) {
				v[i] = min + i*w;
			}
			return v;
		}

		/**
		 * @brief Check if this Axis is empty.
		 */
		[[nodiscard]] 
		bool empty() const noexcept {return bins==0;}

		/**
		 * @brief Get the limits of this Axis.
		 */
		Limit limits() const noexcept {return Limit(min, max);}

		unsigned int bins; 	// The number of equidistant bins. 
		double min;        	// The minimum value spanned by this Axis. 
		double max;        	// The maximum value spanned by this Axis. 
};