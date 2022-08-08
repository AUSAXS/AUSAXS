#pragma once

#include <list>
#include <vector>

#include <data/Atom.h>
#include <hydrate/PlacementStrategy.h>
#include <hydrate/CullingStrategy.h>
#include <hydrate/GridMember.h>
#include <utility/Settings.h>
#include <utility/Exceptions.h>
#include <utility/Axis.h>

// forwards declaration
class Body;

class Grid {
	public:
		/**
		 * @brief Constructor.
		 * 
		 * @param base the base point for the grid.
		 * @param width the distance between each point.
		 * @param bins the number of bins in all dimensions. 
		 */
		Grid(const Axis3D& axes, double width) : Grid(axes, width, setting::grid::ra, setting::grid::rh, setting::grid::placement_strategy, setting::grid::culling_strategy) {}

		/**
		 * @brief Constructor.
		 * 
		 * @param base the base point for the grid.
		 * @param width the distance between each point.
		 * @param bins the number of bins in all dimensions. 
		 * @param radius the radius of each atom.
		 */
		Grid(const Axis3D& axes, double width, int radius) : Grid(axes, width, radius, radius, setting::grid::placement_strategy, setting::grid::culling_strategy) {}

		/**
		 * @brief Constructor.
		 * 
		 * @param base the base point for the grid.
		 * @param width the distance between each point.
		 * @param bins the number of bins in each dimension. 
		 * @param ra the radius of each atom.
		 * @param rh the radius of each water molecule.
		 */
		Grid(const Axis3D& axes, double width, double ra, double rh, setting::grid::PlacementStrategy ps = setting::grid::placement_strategy, setting::grid::CullingStrategy cs = setting::grid::culling_strategy);

		/**
		 * @brief Space-saving constructor. 
		 * 
		 * Construct a new Grid with a volume only slightly larger than that spanned by the input vector. 
		 * 
		 * @param atoms The atoms to be stored in the Grid. 
		 */
		Grid(const std::vector<Atom>& atoms) : Grid(atoms, setting::grid::width, setting::grid::ra, setting::grid::rh, setting::grid::placement_strategy, setting::grid::culling_strategy) {}

		/**
		 * @brief Space-saving constructor. 
		 * 
		 * Construct a new Grid with a volume only slightly larger than that spanned by the input vector. 
		 * 
		 * @param atoms The atoms to be stored in the Grid. 
		 * @param width the distance between each point.
		 * @param bins the number of bins in each dimension. 
		 * @param ra the radius of each atom.
		 * @param rh the radius of each water molecule.
		 */
		Grid(const std::vector<Atom>& atoms, double width, double ra, double rh, setting::grid::PlacementStrategy ps = setting::grid::placement_strategy, setting::grid::CullingStrategy cs = setting::grid::culling_strategy);

		/**
		 * @brief Space-saving constructor. 
		 * 
		 * Construct a new Grid with a volume only slightly larger than that spanned by the input bodies. 
		 * 
		 * @param bodies The bodies to be stored in the Grid. 
		 */
		Grid(const std::vector<Body>& bodies) : Grid(bodies, setting::grid::width, setting::grid::ra, setting::grid::rh, setting::grid::placement_strategy, setting::grid::culling_strategy) {}

		/**
		 * @brief Space-saving constructor. 
		 * 
		 * Construct a new Grid with a volume only slightly larger than that spanned by the input bodies. 
		 * 
		 * @param bodies The bodies to be stored in the Grid. 
		 * @param width the distance between each point.
		 * @param bins the number of bins in each dimension. 
		 * @param ra the radius of each atom.
		 * @param rh the radius of each water molecule.
		 */
		Grid(const std::vector<Body>& bodies, double width, double ra, double rh, setting::grid::PlacementStrategy ps = setting::grid::placement_strategy, setting::grid::CullingStrategy cs = setting::grid::culling_strategy);

		/**
		 * @brief Copy constructor. 
		 */
		Grid(const Grid& grid);

		/** 
		 * @brief Add a set of atoms to the grid. 
		 * @param atoms the set of atoms to add to this grid.
		 */
		template<typename T>
		std::vector<grid::GridMember<T>> add(const std::vector<T>& atoms) {
			static_assert(std::is_base_of<Atom, T>::value, "Argument type must be derivable from Atom!");
			std::vector<grid::GridMember<T>> v(atoms.size());
			size_t index = 0;
			for (const auto& a : atoms) {
				v[index++] = add(a);
			}
			return v;
		}

		/**
		 * @brief Add the contents of a body to this grid.
		 * 
		 * @param body The body to be added. 
		 */
		std::vector<grid::GridMember<Atom>> add(const Body* const body);

		/** 
		 * @brief Add a single atom to the grid. 
		 *        This is a constant-time operation. 
		 */
		grid::GridMember<Atom> add(const Atom& atom, bool expand = false);

		/** 
		 * @brief Add a single hetatom to the grid. 
		 *        This is a constant-time operation. 
		 */
		grid::GridMember<Hetatom> add(const Hetatom& atom, bool expand = false);

		/**
		 * @brief Remove the contents of a body from this grid. 
		 * 
		 * @param body The body to be removed. 
		 */
		void remove(const Body* const body);

		/**
		 * @brief Remove a single atom from the grid.
		 *        This is linear in the number of stored elements. 
		 */
		void remove(const Atom& atom);

		/**
		 * @brief Remove a single water molecule from the grid.
		 *        This is linear in the number of stored elements. 
		 */
		void remove(const Hetatom& atom);

		/**
		 * @brief Remove multiple water molecule from the grid.
		 *        This is linear in the number of stored elements. 
		 */
		void remove(const std::vector<Hetatom>& atom);

		/**
		 * @brief Remove multiple protein atoms from the grid.
		 *        This is linear in the number of stored elements. 
		 */
		void remove(const std::vector<Atom>& atom);

		/**
		 * @brief Remove all waters from the grid.
		 */
		void clear_waters();

		/** 
		 * @brief Expand all member atoms and water molecules into actual spheres based on the radii ra and rh. 
		 */
		void expand_volume();

		/** 
		 * @brief Expand a single member atom into an actual sphere.
		 */
		void expand_volume(grid::GridMember<Atom>& atom);

		/** 
		 * @brief Expand a single member atom into an actual sphere.
		 */
		void expand_volume(grid::GridMember<Hetatom>& atom);

		/** 
		 * @brief Deflate all member atoms and water molecules into actual spheres based on the radii ra and rh. 
		 */
		void deflate_volume();

		/** 
		 * @brief Deflate a single member atom into an actual sphere.
		 */
		void deflate_volume(grid::GridMember<Atom>& atom);

		/** 
		 * @brief Deflate a single member atom into an actual sphere.
		 */
		void deflate_volume(grid::GridMember<Hetatom>& atom);

		/**
		 * @brief Generate a new hydration layer for the grid.
		 * 
		 * @return Pointers to the new water molecules. 
		 */
		std::vector<Hetatom> hydrate();

		/**
		 * @brief Identify possible hydration binding locations for the structure. 
		 * 
		 * @return A list of possible (binx, biny, binz) locations.
		 */
		std::vector<grid::GridMember<Hetatom>> find_free_locs();

		/**
		 * @brief Create the smallest possible box containing the center points of all member atoms.
		 * 
		 * @return Two vectors containing the minimum and maximum coordinates of the box. 
		 */
		std::pair<Vector3<int>, Vector3<int>> bounding_box_index() const;

		static std::pair<Vector3<double>, Vector3<double>> bounding_box(const std::vector<Atom>& atoms);

		/**
		 * @brief Set the radius of all atoms (not water molecules!).
		 * 
		 * @param radius The new radius in Ångström.
		 */
		void set_radius_atoms(double radius);

		/**
		 * @brief Set the radius for all water molecules.
		 * 
		 * @param radius The new radius in Ångström.
		 */
		void set_radius_water(double radius);

		Vector3<int> get_bins() const;

		/**
		 * @brief Get all water molecules from this grid. 
		 */
		std::vector<Hetatom> get_waters() const;

		/**
		 * @brief Get all atoms from this grid. 
		 */
		std::vector<Atom> get_atoms() const;

		/**
		 * @brief Get the total volume spanned by the atoms in this grid in Å^3.
		 *        Water molecules are ignored.  
		 */
		double get_volume();

		/**
		 * @brief Get the width of each bin.
		 */
		double get_width() const {return width;}

		Axis3D get_axes() const {return axes;}

		/**
		 * @brief Convert a vector of bin locations (binx, biny, binz) to a vector of absolute coordinates (x, y, z).
		 * @param v the bin location.
		 * @return The xyz location.
		 */
		template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
		Vector3<double> to_xyz(const Vector3<T>& v) const {
			double x = axes.x.min + width*v[0];
			double y = axes.y.min + width*v[1];
			double z = axes.z.min + width*v[2];
			return {x, y, z};
		}

		/**
		 * @brief Convert a vector of absolute coordinates (x, y, z) to a vector of bin locations.
		 * @param v the xyz location.
		 * @return The bin location. 
		 */
		Vector3<int> to_bins(const Vector3<double>& v) const;

		/**
		 * @brief Get a copy of this Grid. 
		 */
		Grid copy() const;

		/**
		 * @brief Set this Grid equal to another.
		 * 
		 * @param rhs The new Grid. 
		 */
		Grid& operator=(const Grid& rhs);

		/**
		 * @brief Check if this Grid is identical to another. 
		 *        Primarily intended to be used with tests.
		 * 
		 * @param rhs The Grid to compare against. 
		 */
		bool operator==(const Grid& rhs) const;

		std::vector<std::vector<std::vector<char>>> grid; // The actual grid. Datatype is char since we need at least four different values.
		std::list<grid::GridMember<Atom>> a_members; // A list of all member atoms and where they are located.
		std::list<grid::GridMember<Hetatom>> w_members; // A list of all member water molecules and where they are located. 
		unsigned int volume = 0; // The number of bins covered by the members, i.e. the actual volume in the unit (width)^3
		unsigned int ra = 0; // Radius of each atom represented as a number of bins
		unsigned int rh = 0; // Radius of each water molecule represented as a number of bins

	private:
		Axis3D axes;
		double width; // distance between each grid point
		std::unique_ptr<grid::PlacementStrategy> water_placer; // the strategy for placing water molecules
		std::unique_ptr<grid::CullingStrategy> water_culler; // the strategy for culling the placed water molecules

		/** 
		 * @brief Expand a single member atom into an actual sphere.
		 * @param loc The bin location of the atom. 
		 * @param is_water If the atom is a water molecule. Used to determine which marker to use in the grid. 
		 */
		void expand_volume(const Vector3<int>& val, bool is_water);

		/** 
		 * @brief Deflate a single member atom into an actual sphere.
		 * @param loc The bin location of the atom. 
		 * @param is_water If the atom is a water molecule. Used to determine which marker to use in the grid. 
		 */
		void deflate_volume(const Vector3<int>& loc, bool is_water);

		void setup(double width, double ra, double rh, setting::grid::PlacementStrategy ps, setting::grid::CullingStrategy cs);
};