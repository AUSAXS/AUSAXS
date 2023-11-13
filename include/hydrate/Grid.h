#pragma once

#include <hydrate/GridObj.h>
#include <utility/Axis3D.h>
#include <data/DataFwd.h>
#include <io/IOFwd.h>
#include <math/Vector3.h>

#include <list>
#include <vector>
#include <memory>

namespace constants {enum class atom_t;}
namespace grid {
	template<typename T> class GridMember;
	class PlacementStrategy;
	class CullingStrategy;
	class Grid {
		public:
			/**
			 * @brief Constructor.
			 * 
			 * @param axes The axis limits. The width must be set through settings::grid::width.
			 */
			Grid(const Limit3D& axes);

			/**
			 * @brief Space-saving constructor. 
			 * 
			 * Construct a new Grid with a volume only slightly larger than that spanned by the input vector. 
			 * 
			 * @param atoms The atoms to be stored in the Grid. 
			 */
			Grid(const std::vector<data::record::Atom>& atoms);

			/**
			 * @brief Space-saving constructor. 
			 * 
			 * Construct a new Grid with a volume only slightly larger than that spanned by the input bodies. 
			 * 
			 * @param bodies The bodies to be stored in the Grid. 
			 */
			Grid(const std::vector<data::Body>& bodies);

			/**
			 * @brief Copy constructor. 
			 */
			Grid(const Grid& grid);

			/**
			 * @brief Move constructor. 
			 */
			Grid(Grid&& grid) noexcept;

			~Grid();

			/**
			 * @brief Get the atomic radius of an atom in Å.
			 */
			virtual double get_atomic_radius(constants::atom_t atom) const;

			/**
			 * @brief Get the radius of a water molecule in Å.
			 */
			virtual double get_hydration_radius() const;

			/** 
			 * @brief Add a vector of atoms to the grid. 
			 * 		  All added atoms are automatically expanded.
			 * 		  Complexity: O(n) in the number of added atoms.
			 */
			template <typename T, typename = std::enable_if_t<std::is_base_of<data::record::Atom, T>::value>>
			std::vector<grid::GridMember<T>> add(const std::vector<T>& atoms);

			/**
			 * @brief Add the contents of a body to the grid.
			 * 		  All added atoms are automatically expanded.
			 * 		  Complexity: O(n) in the number of added atoms.
			 */
			std::vector<grid::GridMember<data::record::Atom>> add(const data::Body* body);

			/** 
			 * @brief Add a single atom to the grid. 
			 *        Complexity: O(1). 
			 */
			const GridMember<data::record::Atom>& add(const data::record::Atom& atom, bool expand = false);

			/** 
			 * @brief Add a single water to the grid. 
			 *        Complexity: O(1). 
			 */
			const GridMember<data::record::Water>& add(const data::record::Water& atom, bool expand = false);

			/**
			 * @brief Remove the contents of a body from the grid.
			 * 		  All removed atoms are automatically deflated. 
			 * 		  Complexity: O(n) in the number of member atoms.
			 */
			void remove(const data::Body* body);

			/**
			 * @brief Remove atoms as specified by the @a to_remove vector. 
			 * 		  All removed atoms are automatically deflated.
			 * 		  Complexity: O(n) in the number of member atoms.
			 */
			void remove(std::vector<bool>& to_remove);

			/**
			 * @brief Remove a single atom from the grid.
			 *        Complexity: O(n) in the number of member atoms. 
			 */
			void remove(const data::record::Atom& atom);

			/**
			 * @brief Remove a single waters from the grid.
			 *        Complexity: O(n) in the number of member waters. 
			 */
			void remove(const data::record::Water& atom);

			/**
			 * @brief Remove multiple waters from the grid.
			 *        Complexity: O(n) in the number of removed waters. 
			 */
			void remove(const std::vector<data::record::Water>& atom);

			/**
			 * @brief Remove multiple protein atoms from the grid.
			 *        Complexity: O(n) in the number of removed atoms. 
			 */
			void remove(const std::vector<data::record::Atom>& atom);

			/**
			 * @brief Remove all waters from the grid.
			 *        Complexity: O(n) in the number of waters. 
			 */
			void clear_waters();

			/** 
			 * @brief Expand all member atoms and waters into actual spheres based on the radii ra and rh. 
			 * 		  Only expands atoms if they have not already been expanded. 
			 * 		  Complexity: O(n) in the number of unexpanded atoms.
			 */
			void expand_volume();

			/**
			 * @brief Expand all member atoms and waters into actual spheres based on the radii ra and rh.
			 * 		  This method will expand all atoms, regardless of whether they have already been expanded.
			 * 		  Complexity: O(n) in the number of atoms.
			 */
			void force_expand_volume(); 

			/** 
			 * @brief Deflate all member atoms and water molecules into actual spheres based on the radii ra and rh. 
			 * 		  Only deflates atoms if they have been expanded. 
			 * 		  Complexity: O(n) in the number of expanded atoms.
			 */
			void deflate_volume();

			/**
			 * @brief Count the number of atoms in each cluster, and get those with less than \a min atoms.
			 *        This is useful for removing "floating" atoms from e.g. EM map data.
			 * 
			 * @return A vector of booleans indicating whether the atom at the corresponding index is part of a cluster with less than \a min atoms.
			 */
			std::vector<bool> remove_disconnected_atoms(unsigned int min);

			/**
			 * @brief Generate a new hydration layer for the grid.
			 * 		  The layer is generated by the selected strategy.
			 * 		  Complexity: Depends on the strategy. Worst case is O(n) in the number of bins. 
			 * 
			 * @return Pointers to the new water molecules. 
			 */
			std::vector<data::record::Water> hydrate();

			/**
			 * @brief Get the number of bins in each dimension.
			 */
			Vector3<int> get_bins() const;

			/**
			 * @brief Get a copy of all waters in the grid. 
			 * 		  Complexity: O(n) in the number of member waters.
			 */
			std::vector<data::record::Water> get_waters() const;

			/**
			 * @brief Get a copy of all atoms in the grid.
			 * 		  Complexity: O(n) in the number of member atoms.
			 */
			std::vector<data::record::Atom> get_atoms() const;

			/**
			 * @brief Get the total volume spanned by the atoms in this grid in Å^3.
			 * 		  This will trigger the expansion of all unexpanded atoms. 
			 *        Waters do not count towards this volume.
			 * 		  Complexity: O(n) in the number of unexpanded atoms.
			 */
			double get_volume();

			/**
			 * @brief Get the width of each bin.
			 * 		  Complexity: O(1).
			 */
			double get_width() const;

			/**
			 * @brief Get a copy of the axes of the grid.
			 * 		  Complexity: O(1).
			 */
			const Axis3D& get_axes() const {return axes;}

			/**
			 * @brief Create the smallest possible box containing the center points of all member atoms.
			 * 		  Complexity: O(n) in the number of atoms.
			 * 
			 * @return Two vectors containing the minimum and maximum coordinates of the box. 
			 */
			std::pair<Vector3<int>, Vector3<int>> bounding_box_index() const;

			/**
			 * @brief Convert a vector of absolute coordinates (x, y, z) to a vector of bin locations.
			 * 		  Complexity: O(1).
			 */
			Vector3<int> to_bins(const Vector3<double>& v) const;

			/**
			 * @brief Convert a vector of absolute coordinates (x, y, z) to a vector of bin locations.
			 * 		  If the coordinates are outside the grid, they are set to the closest edge.
			 * 		  Complexity: O(1) (but slower than to_bins)
			 */
			Vector3<int> to_bins_bounded(const Vector3<double>& v) const;

			/**
			 * @brief Convert a location in the grid (binx, biny, binz) to a vector of absolute coordinates (x, y, z).
			 * 		  Complexity: O(1).
			 */
			template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
			Vector3<double> to_xyz(const Vector3<T>& v) const {
				return to_xyz(v.x(), v.y(), v.z());
			}

			Vector3<double> to_xyz(int i, int j, int k) const;

			/**
			 * @brief Set this Grid equal to another.
			 * 		  Complexity: O(n) in the number of bins.
			 */
			Grid& operator=(const Grid& rhs);

			/**
			 * @brief Set this Grid equal to another.
			 * 		  Complexity: O(n) in the number of bins.
			 */
			Grid& operator=(Grid&& rhs) noexcept;

			/**
			 * @brief Check if this Grid is identical to another.
			 * 		  Neither the grid contents or members are comparead. 
			 *        Primarily intended to be used with tests.
			 * 		  Complexity: O(1).
			 */
			bool operator==(const Grid& rhs) const;

			/**
			 * @brief Save this Grid as a PDB file.
			 * 		  Complexity: O(n) in the number of bins.
			 */
			void save(const io::File& path) const;

			/**
			 * @brief Convert all bins occupied by atoms to dummy atoms for use in excluded volume calculations.
			 */
			std::vector<Vector3<double>> generate_excluded_volume();

			std::vector<data::record::Atom> get_surface_atoms() const;

			/**
			 * @brief Get the contents of a single bin.
			 */
			const detail::State& index(unsigned int i, unsigned int j, unsigned int k) const;

			/**
			 * @brief Get the center of the grid in bin coordinates.
			 */
			Vector3<int> get_center() const;

			detail::GridObj grid; // The actual grid.
			std::list<GridMember<data::record::Atom>> a_members; // A list of all member atoms and where they are located.
			std::list<GridMember<data::record::Water>> w_members; // A list of all member water molecules and where they are located. 

		protected: // only protected since they are important for testing
			unsigned int volume = 0; // The number of bins covered by the members, i.e. the actual volume in the unit (width)^3

			/** 
			 * @brief Expand a single member atom into an actual sphere.
			 * 		  Only expands atoms if they have not already been expanded. 
			 * 		  Complexity: O(1).
			 */
			void expand_volume(GridMember<data::record::Atom>& atom);

			/** 
			 * @brief Expand a single member water molecule into an actual sphere.
			 * 		  Only expands molecules if they have not already been expanded.
			 * 		  Complexity: O(1).
			 */
			void expand_volume(GridMember<data::record::Water>& atom);

			/** 
			 * @brief Deflate a single member atom into an actual sphere.
			 * 		  Only deflates atoms if they have been expanded.
			 * 		  Complexity: O(1).
			 */
			void deflate_volume(GridMember<data::record::Atom>& atom);

			/** 
			 * @brief Deflate a single member atom into an actual sphere.
			 * 		  Only deflates atoms if they have been expanded.
			 * 		  Complexity: O(1).
			 */
			void deflate_volume(GridMember<data::record::Water>& atom);

			/**
			 * @brief Create the smallest possible box containing all atoms.
			 * 		  Complexity: O(n) in the number of atoms.
			 * 
			 * @return Two vectors containing the minimum and maximum coordinates of the box. 
			 */
			static std::pair<Vector3<double>, Vector3<double>> bounding_box(const std::vector<data::record::Atom>& atoms);

			/**
			 * @brief Identify possible hydration binding locations for the structure. 
			 * 
			 * @return A list of possible (binx, biny, binz) locations.
			 */
			std::vector<GridMember<data::record::Water>> find_free_locs();

		private:
			Axis3D axes;
			std::unique_ptr<PlacementStrategy> water_placer; // the strategy for placing water molecules
			std::unique_ptr<CullingStrategy> water_culler; // the strategy for culling the placed water molecules

			/** 
			 * @brief Expand a single member atom into an actual sphere.
			 * 		  All empty bins within the radius of the atom will be set to either GridObj::A_AREA or GridObj::H_AREA, 
			 * 		  and the center will be set to either GridObj::A_CENTER or GridObj::H_CENTER.
			 * 		  The grid volume is updated accordingly. 
			 * 
			 * @param loc The bin location of the atom. 
			 * @param is_water If the atom is a water molecule. Used to determine which marker to use in the grid. 
			 */
			void expand_volume(const Vector3<int>& loc, bool is_water);

			/** 
			 * @brief Deflate a single member atom from a sphere to a point.
			 * 		  All bins within the radius of the atom will be set to GridObj::EMPTY, and the center
			 * 		  will be set to either GridObj::A_CENTER or GridObj::H_CENTER.
			 * 		  The grid volume is updated accordingly.
			 * 		  
			 * 		  Note that this may cause the grid to be in an inconsistent state without an addtional call to expand_volume, 
			 * 		  since this method does not consider that more than one atom may fill the same bin.
			 * 
			 * @param loc The bin location of the atom. 
			 * @param is_water If the atom is a water molecule. Used to determine which marker to use in the grid. 
			 */
			void deflate_volume(const Vector3<int>& loc, bool is_water);

			void setup();
	};
}