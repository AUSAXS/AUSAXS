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
			StateManager(unsigned int size);

			/**
			 * @brief Mark that the protein atoms of all bodies were internally modified. 
			 */
			void internally_modified_all();

			/**
			 * @brief Mark that the protein atoms of all bodies were externally modified. 
			 */
			void externally_modified_all();

			/**
			 * @brief Mark that the protein atoms of a body was internally modified.
			 * 
			 * @param i index of the body. 
			 */
			void internally_modified(unsigned int i);

			/**
			 * @brief Mark that the protein atoms of a body was externally modified.
			 * 
			 * @param i index of the body. 
			 */
			void externally_modified(unsigned int i);

			/**
			 * @brief Mark that the hydration atoms of a body was modified.
			 * @param i index of the body. 
			 */
			void modified_hydration_layer();

			/**
			 * @brief Reset all marks to false.
			 * ? The awkard name is to avoid accidental collisions with the reset method of smart pointers. 
			 */
			void reset_to_false();

			/**
			 * @brief Get a pointer to the @a ith probe so it can be dispatched to other classes.
			 */
			std::shared_ptr<signaller::Signaller> get_probe(unsigned int i);

			/**
			 * @brief Set the @a ith probe.
			 */
			void set_probe(unsigned int i, std::shared_ptr<signaller::Signaller> probe);

			/**
			 * @brief Get a pointer to all probes.
			 */
			std::vector<std::shared_ptr<signaller::Signaller>> get_probes();

			/**
			 * @brief Get a boolean vector which denotes if the state of a given body was changed. 
			 */
			const std::vector<bool>& get_externally_modified_bodies() const;

			const std::vector<bool>& get_internally_modified_bodies() const;

			/**
			 * @brief Check if a given body has been marked as modified.
			 */
			[[nodiscard]] bool is_externally_modified(unsigned int i);

			/**
			 * @brief Check if a given body has been marked as modified.
			 */
			[[nodiscard]] bool is_internally_modified(unsigned int i);

			/**
			 * @brief Returns true if the hydration layer has been modified, false otherwise. 
			 */
			bool get_modified_hydration() const;

			/**
			 * @brief Get the number of bodies being managed. 
			 */
			[[nodiscard]] std::size_t size() const;

		private:
			std::size_t _size;
			std::vector<bool> _externally_modified;
			std::vector<bool> _internally_modified;
			bool _modified_hydration;
			std::vector<std::shared_ptr<signaller::Signaller>> probes;
	};
}