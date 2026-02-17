// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/HistFwd.h>
#include <data/DataFwd.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <data/symmetry/MoleculeSymmetryFacade.h>
#include <math/MathFwd.h>
#include <io/ExistingFile.h>
#include <utility/observer_ptr.h>
#include <dataset/DatasetFwd.h>
#include <fitter/FitterFwd.h>
#include <grid/GridFwd.h>
#include <hydrate/HydrationFwd.h>
#include <settings/HistogramSettings.h>

#include <string>
#include <vector>
#include <memory>

namespace ausaxs::data {
	/**
	 * @brief A representation of a molecule.
	 * 
	 * A molecule is a collection of bodies. Each body is a collection of atoms.
	 * The hydration layer does not belong to any individual body, and is thus stored in this class. 
	 */
	class Molecule {
		public: 
			Molecule();
			Molecule(const Molecule& other) = delete;
			Molecule& operator=(const Molecule& other) = delete;
			Molecule(Molecule&& other);
			Molecule& operator=(Molecule&& other);
			virtual ~Molecule();

			explicit Molecule(std::vector<Body>&& bodies);
			explicit Molecule(const std::vector<Body>& bodies);

			explicit Molecule(const std::vector<std::string>& input);
			explicit Molecule(const io::File& input);

			/**
			 * @brief Get the distances between each atom.
			 */
			[[nodiscard]] std::unique_ptr<hist::ICompositeDistanceHistogram> get_histogram() const;

			/**
			 * @brief Get the total distance histogram only. 
			 *        This is a slightly faster alternative to get_histogram() when only the total histogram is needed. 
			 */
			[[nodiscard]] std::unique_ptr<hist::DistanceHistogram> get_total_histogram() const;

			/**
			 * @brief Simulate a SAXS dataset based on this molecule.
			 * 
			 * @param add_noise Whether to add noise to the simulated dataset.
			 */
			[[nodiscard]] SimpleDataset simulate_dataset(bool add_noise = true) const;

			/** 
			 * @brief Writes this body to disk.
			 */
			void save(const io::File& path) const;

			/** 
			 * @brief Use an algorithm to generate a new hydration layer for this body. Note that the previous one will be deleted.
			 */
			void generate_new_hydration();

			/**
			 * @brief Calculate the volume of this molecule based on the number of grid bins it spans.
			 * 
			 * @return The volume in Å^3.
			 */
			[[nodiscard]] virtual double get_volume_grid() const;

			/**
			 * @brief Calculate the van der Waals volume of this molecule.
			 * 
			 * @return The volume in Å^3.
			 */
			[[nodiscard]] double get_volume_vdw() const;

			/**
			 * @brief Calculate the excluded volume of this molecule for the currently set excluded volume model. 
			 *
			 * @param d Fitted value of the excluded volume parameter. 
			 *			This must be supplied if the excluded volume has been fitted. 
			 * 
			 * @return The volume in Å^3.
			 */
			[[nodiscard]] double get_volume_exv(double d = 1) const;

			/** 
			 * @brief Calculate the center-mass coordinates.
			 */
			[[nodiscard]] Vector3<double> get_cm(bool include_water = true) const;

			/**
			 * @brief Calculate the atomic molar mass of this molecule in Daltons.
			 * 		  Note that this is just the sum of the molar mass of all atoms.
			 */
			[[nodiscard]] double get_molar_mass() const;

			/**
			 * @brief Get the absolute atomic mass of this entire molecule. 
			 * 		  Note that this is just the sum of the mass of all atoms.  
			 */
			[[nodiscard]] double get_absolute_mass() const;

			/**
			 * @brief Get the excluded volume mass of this molecule.
			 * 		  This is the excluded volume of the molecule times the average protein mass density. 
			 * 
			 * @return The excluded volume mass in Da.
			 */
			[[nodiscard]] double get_excluded_volume_mass() const;

			/**
			 * @brief Get the total atomic charge. 
			 */
			[[nodiscard]] double get_total_atomic_charge() const;

			/**
			 * @brief Get the relative charge density. 
			 */
			[[nodiscard]] double get_relative_charge_density() const;

			/**
			 * @brief Get the relative mass density.
			 */
			[[nodiscard]] double get_relative_mass_density() const;

			/**
			 * @brief Get the relative charge.
			 *        This is the total charge subtracted by the total charge of water of the same volume. 
			 */
			[[nodiscard]] double get_relative_charge() const;

			/**
			 * @brief Get the grid representation. 
			 */
			[[nodiscard]] observer_ptr<grid::Grid> get_grid() const;

			/**
			 * @brief Set the grid representation.
			 */
			void set_grid(grid::Grid&& grid);
			void set_grid(std::unique_ptr<grid::Grid> grid); // @copydoc set_grid(grid::Grid&&)

			/**
			 * @brief Clear the current grid.
			 */
			void clear_grid();

			/**
			 * @brief Remove all waters.
			 */
			void clear_hydration();

			/**
			 * @brief Center this molecule on origo. 
			 */
			void center();

			/**
			 * @brief Get a reference to the body at the given index. 
			 * 		  Complexity: O(1)
			 */
			[[nodiscard]] Body& get_body(unsigned int index);
			[[nodiscard]] const Body& get_body(unsigned int index) const; // @copydoc get_body(unsigned int)

			/**
			 * @brief Get a reference to the bodies of this molecule. 
			 *        Complexity: O(1)
			 */
			[[nodiscard]] std::vector<Body>& get_bodies();
			[[nodiscard]] const std::vector<Body>& get_bodies() const; // @copydoc get_bodies()

			/**
			 * @brief Get a copy of all constituent atoms from the underlying bodies.
			 *        Complexity: O(n)
			 */
			[[nodiscard]] std::vector<data::AtomFF> get_atoms() const;

			/**
			 * @brief Get a copy of all water molecules from the underlying bodies.
			 *        Complexity: O(n)
			 */
			[[nodiscard]] std::vector<data::Water> get_waters() const;

			/**
			 * @brief Get the symmetry facade of this molecule. 
			 *		  This gives access to various symmetry operations, such as getting the full explicit structure.
			 */
			[[nodiscard]] symmetry::detail::MoleculeSymmetryFacade symmetry() const;

			/**
			 * @brief Create a grid and fill it with the atoms of this molecule. 
			 */
			observer_ptr<grid::Grid> create_grid() const;

			/**
			 * @brief Get the radius of gyration of this molecule. 
			 */
			[[nodiscard]] double get_Rg(bool include_waters = true) const;

			/**
			 * @brief Get the number of constituent bodies. 
			 */
			[[nodiscard]] std::size_t size_body() const;

			/**
			 * @brief Get the total number of constituent atoms, excluding hydration. 
			 */
			[[nodiscard]] std::size_t size_atom() const;

			/**
			 * @brief Get the total number of water molecules.
			 */
			[[nodiscard]] std::size_t size_water() const;

			/**
			 * @brief Bind the signaller objects in each body to the histogram manager. 
			 */
			void bind_body_signallers();

			/**
			 * @brief Get the histogram manager of this molecule.
			 */
			[[nodiscard]] observer_ptr<hist::IHistogramManager> get_histogram_manager() const;

			/**
			 * @brief Set the histogram manager of this molecule.
			 */
			void set_histogram_manager(std::unique_ptr<hist::IHistogramManager> manager);
			void set_histogram_manager(settings::hist::HistogramManagerChoice choice);
			void reset_histogram_manager(); // resets to default based on current settings

			/**
			 * @brief Signal that the hydration layer has been modified.
			 */
			void signal_modified_hydration_layer() const;

			/** 
			 * @brief Move the entire molecule by a vector.
			 * @param v the translation vector.
			 */
			void translate(const Vector3<double>& v);

			bool equals_content(const Molecule& other) const;

			/**
			 * @brief Get the hydration generator of this molecule. 
			 */
			[[nodiscard]] observer_ptr<hydrate::HydrationStrategy> get_hydration_generator() const;

			/**
			 * @brief Set the hydration generator of this molecule. 
			 */
			void set_hydration_generator(std::unique_ptr<hydrate::HydrationStrategy> manager);

		private:
			std::vector<Body> bodies;

			// grid is mutable because it is lazily initialized - all methods doing anything but initialization are not const
			mutable std::unique_ptr<grid::Grid> grid; 						// The grid representation of this body
			std::unique_ptr<hist::IHistogramManager> phm;					// The histogram manager of this molecule
			std::unique_ptr<hydrate::HydrationStrategy> hydration_strategy; // The strategy used to generate the hydration layer

			void initialize();
	};
}