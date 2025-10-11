// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/state/DataStateFwd.h>

#include <vector>
#include <memory>

namespace ausaxs::state {
	/**
	 * @brief A state manager which keeps track of changes in each body. 
	 *        This is meant to be used in conjunction with DistanceCalculator, such that it only recalculates what is absolutely necessary. 
	 */
	class StateManager {
		public:
			StateManager(std::size_t size);

			/**
			 * @brief Mark that all atoms of all bodies were internally modified. 
			 */
			void internally_modified_all();

			/**
			 * @brief Mark that all atoms of all bodies were externally modified. 
			 */
			void externally_modified_all();

			/**
			 * @brief Mark that all atoms of the ith body were internally modified.
			 */
			void internally_modified(int i);

			/**
			 * @brief Mark that all atoms of the ith body were externally modified.
			 */
			void externally_modified(int i);

			/**
			 * @brief Mark that the hydration layer was modified.
			 */
			void modified_hydration_layer();

			/**
			 * @brief Mark that the jth symmetry of the ith body was modified.
			 */
			void modified_symmetry(int i, int j);

			/**
			 * @brief Reset all marks to false.
			 * ? The awkard name is to avoid accidental collisions with the reset method of smart pointers. 
			 */
			void reset_to_false();

			/**
			 * @brief Get a pointer to the @a ith probe so it can be dispatched to other classes.
			 */
			std::shared_ptr<signaller::Signaller> get_probe(int i);

			/**
			 * @brief Set the @a ith probe.
			 */
			void set_probe(int i, std::shared_ptr<signaller::Signaller> probe);

			/**
			 * @brief Get a pointer to all probes.
			 */
			std::vector<std::shared_ptr<signaller::Signaller>> get_probes();

			/**
			 * @brief Get a boolean vector which denotes if the external state of a given body was changed. 
			 */
			const std::vector<bool>& get_externally_modified_bodies() const;
			std::vector<bool>& get_externally_modified_bodies(); //< @copydoc get_externally_modified_bodies

			/**
			 * @brief Get a boolean vector which denotes if the internal state of a given body was changed. 
			 */
			const std::vector<bool>& get_internally_modified_bodies() const;
			std::vector<bool>& get_internally_modified_bodies(); //< @copydoc get_internally_modified_bodies

			/**
			 * @brief Get a boolean vector which denotes if the symmetry of a given body was changed. 
			 */
			const std::vector<std::vector<bool>>& get_symmetry_modified_bodies() const;
			std::vector<std::vector<bool>>& get_symmetry_modified_bodies(); //< @copydoc get_symmetry_modified_bodies

			/**
			 * @brief Check if a given body has been marked as modified.
			 */
			[[nodiscard]] bool is_externally_modified(int i) const;

			/**
			 * @brief Check if a given body has been marked as modified.
			 */
			[[nodiscard]] bool is_internally_modified(int i) const;

			/**
			 * @brief Returns true if the jth symmetry of the ith body has been modified, false otherwise. 
			 */
			[[nodiscard]] bool is_modified_symmetry(int i, int j) const;

			/**
			 * @brief Returns true if the hydration layer has been modified, false otherwise. 
			 */
			[[nodiscard]] bool is_modified_hydration() const;

			/**
			 * @brief Returns true if anything has been modified, false otherwise. 
			 */
			[[nodiscard]] bool is_modified() const;

			/**
			 * @brief Get the number of bodies being managed. 
			 */
			[[nodiscard]] std::size_t size() const;

		private:
			std::size_t 					_size;
			std::vector<bool> 				_externally_modified;
			std::vector<bool> 				_internally_modified;
			std::vector<std::vector<bool>> 	_symmetry_modified;
			bool 							_modified_hydration;
			bool 							_modified;
			std::vector<std::shared_ptr<signaller::Signaller>> probes;
	};
}