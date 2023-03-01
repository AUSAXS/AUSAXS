#pragma once

#include <utility/Axis.h>
#include <utility/Limit3D.h>

/**
 * @brief A representation of a 3D coordinate system.
 */
class Axis3D {
	public:
		/**
		 * @brief Default constructor. 
		 */
		Axis3D() noexcept;

		/**
		 * @brief Copy constructor. 
		 */
		Axis3D(const Axis3D& axis) noexcept;

		/**
		 * @brief Constructor.
		 * 
		 * Construct a new Axis3D based on three coordinate axes. 
		 * 
		 * @param x The x-axis. 
		 * @param y The y-axis. 
		 * @param z The z-axis. 
		 */
		Axis3D(const Axis& x, const Axis& y, const Axis& z) noexcept;

		/**
		 * @brief Constructor. 
		 * 
		 * Construct a new Axis3D based on a Limit3D and a bin width.
		 * 
		 * @param limits The limits on each coordinate axis. 
		 * @param width The bin width. 
		 */
		Axis3D(const Limit3D& limits, double width) noexcept;

		/**
		 * @brief Constructor.
		 * 
		 * Construct a new Axis3D based on the minimum and maximum values for each axis, along with the bin width. 
		 */
		Axis3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double width = 1) noexcept;

		/**
		 * @brief Constructor.
		 * 
		 * Construct a new Axis3D based on the minimum and maximum values for each axis, along with the bin width. 
		 * 
		 * @param min A Vector3D containing the minimum values for each coordinate. 
		 * @param max A Vector3D containing the maximum values for each coordinate. 
		 * @param width The bin width. 
		 */
		Axis3D(const Vector3<double>& min, const Vector3<double>& max, double width) noexcept;

		/**
		 * @brief Assignment operator. 
		 * 
		 * Set this object equal to another. 
		 */
		Axis3D& operator=(const Axis3D& rhs) noexcept;

		/**
		 * @brief Equality operator. 
		 * 
		 * Check if this object is equal to another. 
		 */
		bool operator==(const Axis3D& rhs) const noexcept;

		/**
		 * @brief Inequality operator.
		 * 
		 * Check if this object is different from another. 
		 */
		bool operator!=(const Axis3D& rhs) const noexcept;

		/**
		 * @brief Get a string representation of this object. 
		 */
		[[nodiscard]] std::string to_string() const noexcept;

		/**
		 * @brief Stream output operator. 
		 * 
		 * Allows this object to easily be output to a given stream. 
		 */
		friend std::ostream& operator<<(std::ostream& os, const Axis3D& axes) noexcept {os << axes.to_string(); return os;}

		/**
		 * @brief Check if this object is fully initialized. Returns false if any of its Axis are empty.
		 */
		[[nodiscard]] bool empty() const noexcept;

		/**
		 * @brief Recalculate the number of bins for each Axis.
		 * 
		 * @param width The bin width. 
		 */
		void rebin(double width) noexcept;

		Axis x; // The x-axis. 
		Axis y; // The y-axis. 
		Axis z; // The z-axis. 
};