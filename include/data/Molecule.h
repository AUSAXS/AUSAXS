#pragma once

#include <utility/Concepts.h>
#include <utility/view_ptr.h>
#include <hist/HistFwd.h>
#include <io/IOFwd.h>
#include <data/DataFwd.h>

#include <string>
#include <vector>
#include <memory>

class SimpleDataset;
namespace fitter {class Fit;}
namespace grid {class Grid;}
template<numeric T> class Vector3;

namespace data {
	/**
	 * @brief A representation of a molecule.
	 * 
	 * A molecule is a collection of bodies. Each body is a collection of atoms.
	 * The hydration layer does not belong to any individual body, and is thus stored in this class. 
	 */
	class Molecule {
		public: 
			/**
			 * @brief Default constructor. 
			 */
			Molecule() noexcept = default;

			/**
			 * @brief Copy constructor.
			 */
			Molecule(const Molecule& molecule);

			/**
			 * @brief Create a new molecule based on a set of bodies.
			 * 
			 * @param bodies The constituent bodies of this molecule. 
			 */
			explicit Molecule(std::vector<Body>&& bodies);

			/**
			 * @brief Create a new molecule based on a list of input file paths. 
			 * 
			 * @param input A list of paths to the input files. File extensions can be mixed. 
			 */
			explicit Molecule(const std::vector<std::string>& input);

			/**
			 * @brief Create a new molecule based on a single input file path. 
			 * 
			 * @param input Path to the input file. 
			 */
			explicit Molecule(const io::ExistingFile& input);

			/**
			 * @brief Create a new molecule based on a set of bodies.
			 * 
			 * @param bodies The constituent bodies of this molecule. 
			 * @param hydration_atoms The hydration layer. 
			 */
			Molecule(const std::vector<Body>& bodies, const std::vector<record::Water>& hydration_atoms);
			Molecule(const std::vector<Body>& bodies);

			/**
			 * @brief Create a new molecule based on a set of atoms. 
			 * This will only create a single constituent body. 
			 * 
			 * @param molecule_atoms The constituent atoms of this molecule. 
			 * @param hydration_atoms The hydration layer. 
			 */
			Molecule(const std::vector<record::Atom>& molecule_atoms, const std::vector<record::Water>& hydration_atoms);
			Molecule(const std::vector<record::Atom>& molecule_atoms);

			~Molecule();

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
			void save(const io::File& path);

			/** 
			 * @brief Use an algorithm to generate a new hydration layer for this body. Note that the previous one will be deleted.
			 */
			void generate_new_hydration();

			/**
			 * @brief Calculate the volume of this molecule based on its constituent amino acids
			 * 
			 * @return The volume in Å^3.
			 */
			[[nodiscard]] double get_volume_acids() const;

			/**
			 * @brief Calculate the volume of this molecule based on the number of grid bins it spans.
			 * 
			 * @return The volume in Å^3.
			 */
			[[nodiscard]] double get_volume_grid() const;

			/**
			 * @brief Calculate the excluded volume of this molecule based on the number of grid bins it spans.
			 * 		  This is equivalent to get_volume_grid, but with an optional scaling factor applied.
			 * 
			 * @return The volume in Å^3.
			 */
			[[nodiscard]] double get_excluded_volume() const;

			/**
			 * @brief Set the scaling factor for the excluded volume.
			 */
			void set_excluded_volume_scaling(double factor);

			/**
			 * @brief Calculate the volume of this molecule based on the number of C-alpha atoms
			 * 
			 * @return The volume in Å^3.
			 */
			// [[nodiscard]] double get_volume_calpha() const;

			/** 
			 * @brief Calculate the center-mass coordinates.
			 */
			[[nodiscard]] Vector3<double> get_cm() const;

			/**
			 * @brief Calculate the atomic molar mass of this molecule in Daltons.
			 * 		  Note that this is just the sum of the molar mass of all atoms.
			 */
			[[nodiscard]] double molar_mass() const;

			/**
			 * @brief Get the absolute atomic mass of this entire molecule. 
			 * 		  Note that this is just the sum of the mass of all atoms.  
			 */
			[[nodiscard]] double absolute_mass() const;

			/**
			 * @brief Get the total atomic charge. 
			 */
			[[nodiscard]] double total_atomic_charge() const;

			/**
			 * @brief Get the total effective charge. 
			 */
			[[nodiscard]] double total_effective_charge() const;

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
			[[nodiscard]] view_ptr<grid::Grid> get_grid() const;

			/**
			 * @brief Set the grid representation.
			 */
			void set_grid(const grid::Grid&);

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

			/**
			 * @brief Get a reference to the body at the given index. 
			 * 		  Complexity: O(1)
			 */
			[[nodiscard]] const Body& get_body(unsigned int index) const;

			/**
			 * @brief Get a reference to the bodies of this molecule. 
			 *        Complexity: O(1)
			 */
			[[nodiscard]] std::vector<Body>& get_bodies();

			/**
			 * @brief Get a reference to the bodies of this molecule. 
			 *        Complexity: O(1)
			 */
			[[nodiscard]] const std::vector<Body>& get_bodies() const;

			/**
			 * @brief Get a copy of all constituent atoms from the underlying bodies.
			 *        Complexity: O(n)
			 */
			[[nodiscard]] std::vector<record::Atom> get_atoms() const;

			/**
			 * @brief Get a reference to the water molecules of this molecule. 
			 *        Complexity: O(1)
			 */
			[[nodiscard]] std::vector<record::Water>& get_waters();

			/**
			 * @brief Get a reference to the water molecules of this molecule. 
			 *        Complexity: O(1)
			 */
			[[nodiscard]] const std::vector<record::Water>& get_waters() const;

			/**
			 * @brief Get a reference to the specified water molecule.
			 *        Complexity: O(1)
			 */
			[[nodiscard]] record::Water& get_waters(unsigned int i);

			/**
			 * @brief Get a reference to the specified water molecule.
			 *        Complexity: O(1)
			 */
			[[nodiscard]] const record::Water& get_water(unsigned int i) const;

			/**
			 * @brief Create a grid and fill it with the atoms of this molecule. 
			 */
			std::shared_ptr<grid::Grid> create_grid() const;

			/**
			 * @brief Calculate the Debye scattering intensity for this molecule. Does not include hydration atoms. 
			 *        This explicitly calculates each term in the double-sum. For a far more efficient approach, 
			 *        create a ScatteringHistogram and call its equivalent method instead. 
			 */
			[[nodiscard]] std::vector<double> debye_transform() const;

			/**
			 * @brief Get the number of constituent bodies. 
			 */
			[[nodiscard]] unsigned int body_size() const;

			/**
			 * @brief Get the total number of constituent atoms, excluding hydration. 
			 */
			[[nodiscard]] unsigned int atom_size() const;

			/**
			 * @brief Get the total number of water molecules.
			 */
			[[nodiscard]] unsigned int water_size() const;

			/**
			 * @brief Bind the signaller objects in each body to the histogram manager. 
			 */
			void bind_body_signallers();

			/**
			 * @brief Generate a new CRYST1 record for this molecule. 
			 */
			[[deprecated]] void generate_unit_cell();

			/**
			 * @brief Count the number of atoms in each cluster, and remove those with less than \a min atoms.
			 *        This is useful for removing "floating" atoms from e.g. EM map data.
			 */
			[[deprecated]] void remove_disconnected_atoms(double min_percent = 0.05);

			/**
			 * @brief Fit a measurement to this molecule.
			 * 
			 * @param measurement Path to the measurement file to be fitted.
			 */
			[[nodiscard]] std::shared_ptr<fitter::Fit> fit(const io::ExistingFile& measurement) const;

			/**
			 * @brief Get the histogram manager of this molecule.
			 */
			[[nodiscard]] view_ptr<hist::IHistogramManager> get_histogram_manager() const;

			/**
			 * @brief Set the histogram manager of this molecule.
			 */
			void set_histogram_manager(std::unique_ptr<hist::IHistogramManager> manager);

			/**
			 * @brief Signal that the hydration layer has been modified.
			 */
			void signal_modified_hydration_layer() const;

			/**
			 * @brief Subtract the charge of the displaced water molecules from the effective charge of the molecule atoms. 
			 * 
			 * @param scaling The excluded volume scaling factor. Default: 1. 
			 */
			void update_effective_charge(double scaling = 1);

			/** 
			 * @brief Move the entire molecule by a vector.
			 * @param v the translation vector.
			 */
			void translate(const Vector3<double>& v);

			Molecule& operator=(const Molecule& other) = delete;

			bool operator==(const Molecule& other) const;

			bool equals_content(const Molecule& other) const;

		private:
			std::vector<record::Water> hydration_atoms; // Stores the hydration atoms from the generated hydration layer
			std::vector<Body> bodies;           // The constituent bodies

			// the following two variables are only necessary to ensure copying cannot repeat the same work
			bool updated_charge = false;        // True if the effective charge of each atom has been updated to reflect the volume they occupy, false otherwise.
			bool centered = false;              // True if this object has been centered, false otherwise. 

			// grid is mutable because it is lazily initialized - all methods doing anything but initialization are not const
			mutable std::shared_ptr<grid::Grid> grid = nullptr; // The grid representation of this body
			std::shared_ptr<hist::IHistogramManager> phm = nullptr;
			std::shared_ptr<hist::ICompositeDistanceHistogram> histogram = nullptr; // An object representing the distances between atoms
			double excluded_volume_scaling_factor = 1; // The scaling factor for the excluded volume.

			void initialize();
	};
}