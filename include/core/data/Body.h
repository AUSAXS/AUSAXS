// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <data/state/DataStateFwd.h>
#include <data/DataFwd.h>
#include <data/symmetry/SymmetryStorage.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <io/IOFwd.h>
#include <grid/GridFwd.h>
#include <math/MathFwd.h>
#include <hydrate/HydrationFwd.h>

#include <vector>
#include <memory>
#include <optional>

namespace ausaxs::data {
	class Body {
		public:
			Body();
			Body(const Body& body);
			Body(Body&& body) noexcept;
			Body &operator=(const Body& rhs);
			Body &operator=(Body&& rhs) noexcept;
			~Body();

			/** 
			 * @brief Create a new collection of atoms (body) from the input .pdb or .xml file. 
			 * 
			 * @param path path to the input file. 
			 * @param signaller a signalling object to signal changes of state
			 */
			Body(const io::File& path);

			template<AtomVectorFF T>
			explicit Body(T&& atoms);
			explicit Body(const std::vector<Atom>& atoms);

			template<AtomVectorFF T, WaterVector U>
			Body(T&& atoms, U&& waters);

			template<WaterVector T>
			Body(const std::vector<Atom>& atoms, T&& waters);

			/**
			 * @brief Get a reference to the constituent atoms.
			 */
			[[nodiscard]] std::vector<data::AtomFF>& get_atoms();
			[[nodiscard]] const std::vector<data::AtomFF>& get_atoms() const; //< @copydoc get_atoms()

			//? we use std::optional here to explicitly inform the user that it may be empty
			//? principially this is equivalent to returning a potentially zero pointer, but the optional is more explicit
			/**
			 * @brief Get a reference to the hydration atoms.
			 *
			 * @return An optional reference to the hydration atoms. 
			 *   	   This may be empty if the hydration model does not support explicit dummy atoms. 
			 */
			[[nodiscard]] std::optional<std::reference_wrapper<std::vector<data::Water>>> get_waters();
			[[nodiscard]] std::optional<std::reference_wrapper<const std::vector<data::Water>>> get_waters() const; //< @copydoc get_waters()

			/**
			 * @brief Set the hydration object.
			 */
			void set_hydration(std::unique_ptr<hydrate::Hydration> hydration);

			/**
			 * @brief Clear the hydration layer.
			 */
			void clear_hydration();
	
			/**
			 * @brief Get a reference to the specified atom.
			 */
			[[nodiscard]] data::AtomFF& get_atom(unsigned int index);
			[[nodiscard]] const data::AtomFF& get_atom(unsigned int index) const;

			/** 
			 * @brief Calculate the center-of-mass coordinates for the body.
			 */
			[[nodiscard]] Vector3<double> get_cm() const;
			
			/**
			 * @brief Calculate the van der Waals volume of this body.
			 * 
			 * @return The volume in Å^3.
			 */
			[[nodiscard]] double get_volume_vdw() const;

			/**
			 * @brief Calculate the molar mass of this body in Daltons.
			 */
			[[nodiscard]] double get_molar_mass() const;

			/**
			 * @brief Get the absolute mass of this body in kg.
			 */
			[[nodiscard]] double get_absolute_mass() const;

			/**
			 * @brief Get the total atomic charge of this body.
			 */
			[[nodiscard]] double get_total_atomic_charge() const;

			/** 
			 * @brief Translate all atoms by a given vector.
			 */
			void translate(Vector3<double> v);

			/**
			 * @brief Rotate all atoms by a given rotation matrix.
			 */
			void rotate(const Matrix<double>& R);

			/**
			 * @brief Check if this object is equal to another based on their unique ID. 
			 */
			bool operator==(const Body& rhs) const;

			/**
			 * @brief Check if the content of this object is equal to another, disregarding their unique ID. 
			 */
			bool equals_content(const Body& rhs) const;

			/**
			 * @brief Get the unique identifier of this Body. 
			 */
			[[nodiscard]] int get_uid() const;

			/**
			 * @brief Get the total number of constituent atoms, excluding waters. 
			 */
			[[nodiscard]] std::size_t size_atom() const;

			/**
			 * @brief Get the total number of water molecules.
			 */
			[[nodiscard]] std::size_t size_water() const;
			
			/**
			 * @brief Get the number of symmetries defined for this body.
			 */
			[[nodiscard]] std::size_t size_symmetry() const;

			/**
			 * @brief Get the number of symmetry duplicates of this body.
			 * 		  This accounts for repeating symmetries. 
			 */
			[[nodiscard]] std::size_t size_symmetry_total() const;

			/**
			 * @brief Access the symmetry operations of this body.
			 */
			[[nodiscard]] symmetry::detail::BodySymmetryFacade<Body> symmetry();
			[[nodiscard]] symmetry::detail::BodySymmetryFacade<const Body> symmetry() const; //< @copydoc symmetry()

			/**
			 * @brief Register a probe (listener) to this object, which will be notified of state changes. 
			 */
			void register_probe(std::shared_ptr<signaller::Signaller> signal);

			/**
			 * @brief Get the signaller object for this body. 
			 */
			std::shared_ptr<signaller::Signaller> get_signaller() const;

		private:
			std::vector<data::AtomFF> 					atoms;
			std::unique_ptr<hydrate::Hydration> 		hydration;
			std::unique_ptr<symmetry::SymmetryStorage> 	symmetries;

			int uid;
			inline static unsigned int uid_counter = 0;

			// The signalling object to signal a change of state. The default doesn't do anything, and must be overriden by a proper Signaller object.  
			std::shared_ptr<signaller::Signaller> signal;

			void initialize();

		friend class symmetry::detail::BodySymmetryFacade<Body>;
		friend class symmetry::detail::BodySymmetryFacade<const Body>;
	};
	static_assert(std::is_move_constructible_v<Body>);
}