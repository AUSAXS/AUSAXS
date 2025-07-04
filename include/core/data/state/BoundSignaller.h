// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/state/DataStateFwd.h>
#include <data/state/Signaller.h>

namespace ausaxs::signaller {
    /**
     * @brief A small probe for signalling changes which can be dispatched to other classes. 
     */
    class BoundSignaller : public Signaller {
        public: 
			BoundSignaller(const BoundSignaller& rhs) = default;
			BoundSignaller(BoundSignaller&& rhs) noexcept = default;
			~BoundSignaller() override = default;

            BoundSignaller(unsigned int id, state::StateManager* const owner);

            void modified_external() const override;
            void modified_internal() const override;
            void modified_symmetry(int i) const override;
            void modified_hydration() const override;
            void set_symmetry_size(std::size_t size) const override;

            /**
             * @brief Get the id of this signaller. 
             */
            unsigned int get_id() const;

        private: 
            state::StateManager* const owner;
            unsigned int id;
    };
}