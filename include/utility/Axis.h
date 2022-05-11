	#pragma once

	#include <vector>
	#include <array>
	#include <initializer_list>
	#include <math.h>

	/**
	 * @brief \class Limit
	 * 
	 * A simple representation of a limited span of values. 
	 */
	class Limit {
		public:
			/**
			 * @brief Default constructor.
			 */
			Limit() : min(0), max(0) {}

			/**
			 * @brief Constructor. 
			 * 
			 * @param min Minimum value of the limit. 
			 * @param max Maximum value of the limit. 
			 */
			Limit(double min, double max) : min(min), max(max) {}

			/**
			 * @brief Get the interval spanned by this limit. 
			 */
			double span() const {return max-min;}

			/**
			 * @brief Equality operator. Check if this Limit is equal to another.
			 */
			bool operator==(const Limit& rhs) const {return min == rhs.min && max == rhs.max;}

			/**
			 * @brief Inequality operator.
			 * 
			 * Check if this object is different from another. 
			 */
			bool operator!=(const Limit& rhs) const {return !operator==(rhs);}

			/**
			 * @brief Stream output operator. 
			 * 
			 * Allows this object to easily be output to a given stream. 
			 */
			friend std::ostream& operator<<(std::ostream& os, const Limit& axis) {os << axis.to_string(); return os;}

			/**
			 * @brief Get a string representation of this object. 
			 */
			std::string to_string() const {return "Limits: (" + std::to_string(min) + ", " + std::to_string(max) + ")";}

			/**
			 * @brief Check if this Limit is empty (min == max == 0).
			 */
			[[nodiscard]] 
			bool empty() const {return min == 0 && max == 0;}

			double min; // The minimum value of this limit. 
			double max; // The maximum value of this limit. 
	};

	/**
	 * @brief \class Limit3D
	 * 
	 * A representation of a 3-dimensional limited span of values. 
	 */
	class Limit3D {
		public:
			Limit3D() {}
			/**
			 * @brief Constructor.
			 * 
			 * Construct a new 3D limit based on three 1D limits. 
			 * 
			 * @param x The limit on the x-axis. 
			 * @param y The limit on the y-axis. 
			 * @param z The limit on the z-axis. 
			 */
			Limit3D(const Limit& x, const Limit& y, const Limit& z) : x(x), y(y), z(z) {}

			/**
			 * @brief Constructor. 
			 * 
			 * Construct a new 3D limit based on the minimum and maximum values along each axis. 
			 */
			Limit3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) : x(xmin, xmax), y(ymin, ymax), z(zmin, zmax) {}

			/**
			 * @brief Check if this object is fully initialized. Returns false if any of its Limits are empty.
			 */
			[[nodiscard]] 
			bool empty() const {return x.empty() || y.empty() || z.empty();}

			Limit x; // The 1D limit on the x-axis. 
			Limit y; // The 1D limit on the y-axis. 
			Limit z; // The 1D limit on the z-axis. 
	};

	/**
	 * @brief \class Axis
	 * 
	 * A representation of an axis. 
	 */
	class Axis {
		public:
			/**
			 * @brief Default constructor. 
			 */
			Axis() : bins(0), min(0), max(0) {}

			/**
			 * @brief Constructor. 
			 * 
			 * Construct a new axis based on a \class Limit and the number of bins. 
			 * 
			 * @param limits The limits on the axis. 
			 * @param bins The number of equidistant bins. 
			 */
			Axis(const Limit& limits, int bins) : bins(bins), min(limits.min), max(limits.max) {}

			/**
			 * @brief Constructor. 
			 * 
			 * Construct a new axis based on explicit minimum and maximum values, along with the number of bins. 
			 * 
			 * @param bins The number of equidistant bins. 
			 * @param xmin The minimum value spanned by this axis. 
			 * @param xmax The maximum value spanned by this axis. 
			 */
			Axis(const int bins, const double xmin, const double xmax) : bins(bins), min(xmin), max(xmax)  {}

			/**
			 * @brief List initializer. 
			 * 
			 * Initialize this Axis with an \class initializer_list. 
			 */
			Axis& operator=(std::initializer_list<double> list) {
				std::vector<double> d = list;
				bins = std::round(d[0]); 
				min = d[1];
				max = d[2];
				return *this;
			}

			/**
			 * @brief Get a string representation of this object. 
			 */
			std::string to_string() const {
				return "Axis: (" + std::to_string(min) + ", " + std::to_string(max) + ") with " + std::to_string(bins) + " bins";
			}

			/**
			 * @brief Equality operator. 
			 * 
			 * Check if this object is equal to another. 
			 */
			bool operator==(const Axis& rhs) const {
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
			bool operator!=(const Axis& rhs) const {return !operator==(rhs);}

			/**
			 * @brief Stream output operator. 
			 * 
			 * Allows this object to easily be output to a given stream. 
			 */
			friend std::ostream& operator<<(std::ostream& os, const Axis& axis) {os << axis.to_string(); return os;}

			/**
			 * @brief Get the bin width. 
			 */
			double width() const {return (max-min)/bins;}

			/**
			 * @brief Get the bin width.
			 */
			double step() const {return width();}

			/**
			 * @brief Check if this Axis is empty.
			 */
			[[nodiscard]] 
			bool empty() const {return bins==0;}

			unsigned int bins; // The number of equidistant bins. 
			double min;        // The minimum value spanned by this Axis. 
			double max;        // The maximum value spanned by this Axis. 
	};

	/**
	 * @brief \class Axis3D
	 * 
	 * A representation of a 3D coordinate system.
	 */
	class Axis3D {
		public:
			/**
			 * @brief Default constructor. 
			 */
			Axis3D() {}

			/**
			 * @brief Copy constructor. 
			 */
			Axis3D(const Axis3D& axis) : x(axis.x), y(axis.y), z(axis.z) {}

			/**
			 * @brief Constructor.
			 * 
			 * Construct a new Axis3D based on three coordinate axes. 
			 * 
			 * @param x The x-axis. 
			 * @param y The y-axis. 
			 * @param z The z-axis. 
			 */
			Axis3D(const Axis& x, const Axis& y, const Axis& z) : x(x), y(y), z(z) {}

			/**
			 * @brief Constructor. 
			 * 
			 * Construct a new Axis3D based on a \class Limit3D and a bin width.
			 * 
			 * @param limits The limits on each coordinate axis. 
			 * @param width The bin width. 
			 */
			Axis3D(const Limit3D& limits, double width) : x(limits.x, limits.x.span()/width), y(limits.y, limits.y.span()/width), z(limits.z, limits.z.span()/width) {}

			/**
			 * @brief Constructor.
			 * 
			 * Construct a new Axis3D based on the minimum and maximum values for each axis, along with the number of equidistant bins. 
			 */
			Axis3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int bins) : x(bins, xmin, xmax), y(bins, ymin, ymax), z(bins, zmin, zmax) {}

			/**
			 * @brief Constructor.
			 * 
			 * Construct a new Axis3D based on the minimum and maximum values for each axis, along with the bin width. 
			 */
			Axis3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double width) : x((xmax-xmin)/width, xmin, xmax), y((ymax-ymin)/width, ymin, ymax), z((zmax-zmin)/width, zmin, zmax) {}

			/**
			 * @brief Constructor.
			 * 
			 * Construct a new Axis3D based on the minimum and maximum values for each axis, along with the bin width. 
			 * 
			 * @param min A 3D vector containing the minimum values for each coordinate. 
			 * @param max A 3D vector containing the maximum values for each coordinate. 
			 * @param width The bin width. 
			 */
			Axis3D(const std::vector<int>& min, const std::vector<int>& max, double width) : x((max[0]-min[0])/width, min[0], max[0]), y((max[1]-min[1])/width, min[1], max[1]), z((max[2]-min[2])/width, min[2], max[2]) {}

			/**
			 * @brief Assignment operator. 
			 * 
			 * Set this object equal to another. 
			 */
			Axis3D& operator=(const Axis3D& rhs) {
			x = rhs.x;
			y = rhs.y;
			z = rhs.z;
			return *this;
			}

			/**
			 * @brief Equality operator. 
			 * 
			 * Check if this object is equal to another. 
			 */
			bool operator==(const Axis3D& rhs) const {
			if (x != rhs.x) {return false;}
			if (y != rhs.y) {return false;}
			if (z != rhs.z) {return false;}
			return true;
			}

			/**
			 * @brief Inequality operator.
			 * 
			 * Check if this object is different from another. 
			 */
			bool operator!=(const Axis3D& rhs) const {return !operator==(rhs);}

			/**
			 * @brief Get a string representation of this object. 
			 */
			std::string to_string() const {
			return "Axes: \n\tx: " + x.to_string() + "\n\ty: " + y.to_string() + "\n\tz: " + z.to_string(); 
			}

			/**
			 * @brief Stream output operator. 
			 * 
			 * Allows this object to easily be output to a given stream. 
			 */
			friend std::ostream& operator<<(std::ostream& os, const Axis3D& axes) {os << axes.to_string(); return os;}

			/**
			 * @brief Check if this object is fully initialized. Returns false if any of its Axis are empty.
			 */
			[[nodiscard]] 
			bool empty() const {return x.empty() || y.empty() || z.empty();}

			Axis x; // The x-axis. 
			Axis y; // The y-axis. 
			Axis z; // The z-axis. 
	};