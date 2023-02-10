#pragma once

#include <string>
#include <vector>
#include <memory>

#include <data/Body.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <hist/detail/HistogramManagerFactory.h>
#include <utility/SimpleDataset.h>
#include <fitter/Fit.h>

class Protein {
	public: 
		/**
		 * @brief Default constructor. 
		 */
		Protein() noexcept {}

		/**
		 * @brief Copy constructor.
		 * 
		 * Proteins cannot be copied. 
		 */
		Protein(const Protein& protein) = delete;

		/**
		 * @brief Move constructor.
		 */
		Protein(Protein&& protein) noexcept;

		/**
		 * @brief Constructor.
		 * 
		 * Create a new protein based on a set of bodies.
		 * 
		 * @param bodies The constituent bodies of this protein. 
		 * @param hydration_atoms The hydration layer. 
		 */
		explicit Protein(const std::vector<Body>& bodies, const std::vector<Water>& hydration_atoms = {});

		/**
		 * @brief Constructor.
		 * 
		 * Create a new protein based on a set of atoms. 
		 * This will only create a single constituent body. 
		 * 
		 * @param protein_atoms The constituent atoms of this protein. 
		 * @param hydration_atoms The hydration layer. 
		 */
		explicit Protein(const std::vector<Atom>& protein_atoms, const std::vector<Water>& hydration_atoms = {});

		/**
		 * @brief Constructor. 
		 * 
		 * Create a new protein based on a set of atom vectors. Each vector defines a constituent body. 
		 * 
		 * @param protein_atoms The constituent atoms of each body. 
		 * @param hydration_atoms The hydration layer. 
		 */
		explicit Protein(const std::vector<std::vector<Atom>>& protein_atoms, const std::vector<Water>& hydration_atoms = {});

		/**
		 * @brief Constructor. 
		 * 
		 * Create a new protein based on a list of input file paths. 
		 * 
		 * @param input A list of paths to the input files. File extensions can be mixed. 
		 */
		explicit Protein(const std::vector<std::string>& input);

		/**
		 * @brief Constructor.
		 * 
		 * Create a new protein based on a single input file path. 
		 * 
		 * @param input Path to the input file. 
		 */
		explicit Protein(std::string input);

		/**
		 * @brief Get the distances between each atom.
		 */
		hist::ScatteringHistogram get_histogram();

		/**
		 * @brief Get the total distance histogram only. 
		 *        This is a slightly faster alternative to get_histogram() when only the total histogram is needed. 
		 */
		hist::Histogram get_total_histogram();

		/**
		 * @brief Simulate a SAXS dataset based on this protein.
		 * 
		 * @param add_noise Whether to add noise to the simulated dataset.
		 */
		SimpleDataset simulate_dataset(bool add_noise = true);

		/** 
		 * @brief Writes this body to disk.
		 */
		void save(std::string path);

		/** 
		 * @brief Use an algorithm to generate a new hydration layer for this body. Note that the previous one will be deleted.
		 */
		void generate_new_hydration();

		/**
		 * @brief Calculate the volume of this protein based on its constituent amino acids
		 * 
		 * @return The volume in Å^3.
		 */
		double get_volume_acids() const;

		/**
		 * @brief Calculate the volume of this protein based on the number of grid bins it spans.
		 * 
		 * @return The volume in Å^3.
		 */
		double get_volume_grid();

		/**
		 * @brief Calculate the volume of this protein based on the number of C-alpha atoms
		 * 
		 * @return The volume in Å^3.
		 */
		double get_volume_calpha() const;

		/** 
		 * @brief Calculate the center-mass coordinates.
		 */
		Vector3<double> get_cm() const;

		/**
		 * @brief Calculate the molar mass of this protein in Daltons.
		 */
		[[nodiscard]] double molar_mass() const;

		/**
		 * @brief Get the absolute mass.
		 * 
		 * @return The mass in kg.
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
		double get_relative_charge_density();

		/**
		 * @brief Get the relative mass density.
		 */
		double get_relative_mass_density();

		/**
		 * @brief Get the relative charge.
		 *        This is the total charge subtracted by the total charge of water of the same volume. 
		 */
		double get_relative_charge();

		/**
		 * @brief Get the grid representation. 
		 */
		std::shared_ptr<Grid> get_grid();

		/**
		 * @brief Set the grid representation. 
		 */
		void set_grid(const Grid& grid);

		/**
		 * @brief Clear the current grid.
		 */
		void clear_grid();

		/**
		 * @brief Remove all waters.
		 */
		void clear_hydration();

		/**
		 * @brief Create a binding point between two bodies.
		 *        This binding will be a constraint for rigid-body optimization. 
		 */
		void bind();

		/**
		 * @brief Center this protein on origo. 
		 */
		void center();

		/**
		 * @brief Get a copy of all constituent atoms from the underlying bodies.
		 *        Complexity: O(n)
		 */
		std::vector<Atom> atoms() const;

		/**
		 * @brief Get a reference to the water molecules of this protein. 
		 *        Complexity: O(1)
		 */
		std::vector<Water>& waters();

		/**
		 * @brief Get a reference to the water molecules of this protein. 
		 *        Complexity: O(1)
		 */
		const std::vector<Water>& waters() const;

		/**
		 * @brief Create a grid and fill it with the atoms of this protein. 
		 */
		std::shared_ptr<Grid> create_grid();

		/**
		 * @brief Calculate the Debye scattering intensity for this protein. Does not include hydration atoms. 
		 *        This explicitly calculates each term in the double-sum. For a far more efficient approach, 
		 *        create a ScatteringHistogram and call its equivalent method instead. 
		 */
		std::vector<double> calc_debye_scattering_intensity();

		/**
		 * @brief Get the number of constituent bodies. 
		 */
		size_t body_size() const;

		/**
		 * @brief Get the total number of constituent atoms, excluding hydration. 
		 */
		size_t atom_size() const;

		/**
		 * @brief Get the total number of constituent atoms, excluding hydration. Equivalent to \a body_size. 
		 */    
		size_t size() const {return atom_size();}

		/**
		 * @brief Bind the signaller objects in each body to the histogram manager. 
		 */
		void bind_body_signallers();

		/**
		 * @brief Generate a new CRYST1 record for this protein. 
		 */
		void generate_unit_cell();

		/**
		 * @brief Count the number of atoms in each cluster, and remove those with less than \a min atoms.
		 *        This is useful for removing "floating" atoms from e.g. EM map data.
		 */
		void remove_disconnected_atoms(unsigned int min = 10);

		/**
		 * @brief Fit a measurement to this protein.
		 * 
		 * @param measurement Path to the measurement file to be fitted.
		 */
		std::shared_ptr<Fit> fit(std::string measurement);

		std::shared_ptr<hist::HistogramManager> get_histogram_manager() const;

		/**
		 * @brief Create a new histogram manager of type \a T.
		 * 
		 * @tparam T: Manager to create. Only the PartialHistogramManagerMT is intended for production. Options:
		 *         - HistogramManager: A simple manager that recalculates the entire histogram every time.
		 *         - HistogramManagerMT: A multithreaded implementation of the simple manager.
		 *         - PartialHistogramManager: A smart manager that only recalculates the parts of the histogram that are needed.
		 *         - PartialHistogramManagerMT: A multithreaded implementation of the partial manager.
		 */
		template<hist::detail::HistogramManagerType T>
		void set_histogram_manager() {phm = hist::HistogramManagerFactory::create<T>(this);}

		std::vector<Water> hydration_atoms; // Stores the hydration atoms from the generated hydration layer
		std::vector<Body> bodies;           // The constituent bodies
		bool updated_charge = false;        // True if the effective charge of each atom has been updated to reflect the volume they occupy, false otherwise
		bool centered = false;              // True if this object is centered, false otherwise. 
	private:
		std::shared_ptr<Grid> grid = nullptr; // The grid representation of this body
		std::shared_ptr<hist::HistogramManager> phm = nullptr;
		std::shared_ptr<hist::ScatteringHistogram> histogram = nullptr; // An object representing the distances between atoms

		/** 
		 * @brief Move the entire protein by a vector.
		 * @param v the translation vector.
		 */
		void translate(const Vector3<double>& v);

		/**
		 * @brief Subtract the charge of the displaced water molecules from the effective charge of the protein atoms. 
		 */
		void update_effective_charge();
};