// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <grid/detail/GridObj.h>
#include <grid/detail/GridInternalFwd.h>
#include <grid/detail/GridExcludedVolume.h>
#include <utility/Axis3D.h>
#include <utility/TypeTraits.h>
#include <data/DataFwd.h>
#include <io/IOFwd.h>
#include <math/Vector3.h>
#include <constants/ConstantsFwd.h>
#include <form_factor/FormFactorType.h>

#include <vector>
#include <span>

namespace ausaxs::grid {
	class Grid {
		struct private_ctr {explicit private_ctr() = default;};
		public:
			/**
			 * @brief Initialize a new grid of the given size with the given cell width. 
			 * 		  This can only be used internally by the Grid class. 
			 */
			Grid(const Axis3D& axes, private_ctr);

			/**
			 * @brief Initialize a new grid of the given size. The cell width is controlled by the settings::grid::cell_width variable.
			 */
			Grid(const Limit3D& axes);

			/**
			 * @brief Space-saving constructor. 
			 * 
			 * Construct a new Grid with a volume only slightly larger than that spanned by the input vector. 
			 * 
			 * @param atoms The atoms to be stored in the Grid. 
			 */
			Grid(const std::vector<data::AtomFF>& atoms);

			/**
			 * @brief Space-saving constructor. 
			 * 
			 * Construct a new Grid with a volume only slightly larger than that spanned by the input bodies. 
			 * 
			 * @param bodies The bodies to be stored in the Grid. 
			 */
			Grid(const std::vector<data::Body>& bodies);

			Grid(const Grid& grid);
			Grid(Grid&& grid) noexcept;
			virtual ~Grid();

			/**
			 * @brief Create a new Grid from a reference file.
			 * 
			 * @param path The path to the reference file. 
			 */
			[[nodiscard]] static std::unique_ptr<Grid> create_from_reference(const io::ExistingFile& path, const data::Molecule& molecule);

			/**
			 * @brief Get the atomic radius of an atom in Å.
			 */
			[[nodiscard]] virtual double get_atomic_radius(form_factor::form_factor_t atom) const;

			/**
			 * @brief Get the radius of a water molecule in Å.
			 */
			[[nodiscard]] virtual double get_hydration_radius() const;

			/**
			 * @brief Add the contents of a body to the grid.
			 *		  Waters are ignored and must be added separately.
			 * 		  Complexity: O(n) in the number of added atoms.
			 */
			std::span<grid::GridMember<data::AtomFF>> add(const data::Body& body, bool expand = true);

			/**
			 * @brief Add waters to the grid.
			 */
			std::span<grid::GridMember<data::Water>> add(const std::vector<data::Water>& waters, bool expand = true);

			/**
			 * @brief Add a single water to the grid.
			 * 		  Complexity: O(1).
			 */
			grid::GridMember<data::Water>& add(const data::Water& water, bool expand = true);

			/**
			 * @brief Remove the contents of a body from the grid.
			 * 		  All removed atoms are automatically deflated. 
			 * 		  Complexity: O(n) in the number of member atoms.
			 */
			void remove(const data::Body& body);

			/**
			 * @brief Remove the given waters. 
			 * 		  All removed waters are automatically deflated.
			 * 		  Complexity: O(n) in the number of waters.
			 */
			void remove_waters(const std::vector<bool>& to_remove);

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
			 * @brief Get the number of bins in each dimension.
			 */
			[[nodiscard]] Vector3<int> get_bins() const;

			/**
			 * @brief Get the total volume spanned by the atoms in this grid in Å^3.
			 * 		  This will trigger the expansion of all unexpanded atoms. 
			 *        Waters may count towards this volume dependending on the settings.
			 * 		  Complexity: O(n) in the number of unexpanded atoms.
			 */
			[[nodiscard]] double get_volume();

			/**
			 * @brief Get the total volume spanned by the atoms in this grid as a number of bins.
			 * 		  This will trigger the expansion of all unexpanded atoms. 
			 *        Waters may count towards this volume dependending on the settings.
			 * 		  Complexity: O(n) in the number of unexpanded atoms.
			 */
			[[nodiscard]] int get_volume_bins() const {return volume;}

			/**
			 * @brief Get the width of each bin.
			 * 		  Complexity: O(1).
			 */
			[[nodiscard]] double get_width() const;

			/**
			 * @brief Get a copy of the axes of the grid.
			 * 		  Complexity: O(1).
			 */
			[[nodiscard]] const Axis3D& get_axes() const {return axes;}

			/**
			 * @brief Convert a vector of absolute coordinates (x, y, z) to a vector of bin locations.
			 * 		  Complexity: O(1).
			 */
			[[nodiscard]] Vector3<int> to_bins(const Vector3<double>& v) const;

			/**
			 * @brief Convert a vector of absolute coordinates (x, y, z) to a vector of bin locations.
			 * 		  If the coordinates are outside the grid, they are set to the closest edge.
			 * 		  Complexity: O(1) (but slower than to_bins)
			 */
			[[nodiscard]] Vector3<int> to_bins_bounded(const Vector3<double>& v) const;

			/**
			 * @brief Convert a location in the grid (binx, biny, binz) to a vector of absolute coordinates (x, y, z).
			 * 		  Complexity: O(1).
			 */
			template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
			[[nodiscard]] Vector3<double> to_xyz(const Vector3<T>& v) const {
				return to_xyz(v.x(), v.y(), v.z());
			}

			[[nodiscard]] Vector3<double> to_xyz(int i, int j, int k) const; //< @copydoc to_xyz(const Vector3<T>& v)

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
			 * 		  This will expand all atoms in the grid.
			 */
			[[nodiscard]] virtual exv::GridExcludedVolume generate_excluded_volume();

			[[nodiscard]] std::vector<data::AtomFF> get_surface_atoms() const;

			/**
			 * @brief Get the contents of a single bin.
			 */
			[[nodiscard]] const detail::State& index(int i, int j, int k) const;

			/**
			 * @brief Get the center of the grid in bin coordinates.
			 */
			[[nodiscard]] Vector3<int> get_center() const;

			/**
			 * @brief Get the atoms in the grid.
			 *		  Complexity: O(N).
			 */
			[[nodiscard]] std::vector<data::AtomFF> get_atoms();

			/**
			 * @brief Get the water molecules in the grid.
			 *		  Complexity: O(N).
			 */
			[[nodiscard]] std::vector<data::Water> get_waters();

			/**
			 * @brief Add a value to the volume of the grid.
			 */
			void add_volume(int value);

			/**
			 * @brief Convert a x bin index to a real x coordinate.
			 * 		  Complexity: O(1).
			 */
			[[nodiscard]] double to_x(int i) const;

			/**
			 * @brief Convert a y bin index to a real y coordinate.
			 * 		  Complexity: O(1).
			 */
			[[nodiscard]] double to_y(int j) const;
 
			/**
			 * @brief Convert a z bin index to a real z coordinate.
			 * 		  Complexity: O(1).
			 */
			[[nodiscard]] double to_z(int k) const;

			/**
			 * @brief Create the smallest possible box containing the center points of all member atoms.
			 * 		  Complexity: O(n) in the number of atoms.
			 * 
			 * @return Two vectors containing the minimum and maximum coordinates of the box. 
			 */
			[[nodiscard]] std::pair<Vector3<int>, Vector3<int>> bounding_box_index(bool include_waters = false) const;

			/**
			 * @brief Create the smallest possible box containing all atoms.
			 * 		  Complexity: O(n) in the number of atoms.
			 * 
			 * @return Two vectors containing the minimum and maximum coordinates of the box. 
			 */
			[[nodiscard]] static std::pair<Vector3<double>, Vector3<double>> 
				bounding_box(const std::vector<data::AtomFF>& atoms);

			// @copydoc bounding_box(const std::vector<data::AtomFF>& atoms)
			[[nodiscard]] static std::pair<Vector3<double>, Vector3<double>> 
				bounding_box(const std::vector<data::Water>& atoms);

			detail::GridObj grid; // The actual grid.
			std::vector<GridMember<data::AtomFF>> a_members; // The member atoms and where they are located.
			std::vector<GridMember<data::Water>>  w_members; // The member water molecules and where they are located. 
			std::unordered_map<int, int> body_start; 		 // The starting index of each body in the a_members vector. 

		protected:
			int volume = 0; // The number of bins covered by the members, i.e. the actual volume in the unit (width)^3
				
		private:
			Axis3D axes;

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

			/**
			 * @brief Remove atoms as specified by the @a to_remove vector. 
			 * 		  All removed atoms are automatically deflated.
			 * 		  Complexity: O(n) in the number of member atoms.
			 */
			void remove(const std::vector<bool>& to_remove);

			void setup();
	};
	static_assert(supports_nothrow_move_v<Grid>, "Grid should be nothrow move constructible");
}