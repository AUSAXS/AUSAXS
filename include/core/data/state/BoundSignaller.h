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

            virtual void external_change() const override;

            virtual void internal_change() const override;

            virtual void symmetry_changed() const override;

            /**
             * @brief Get the id of this signaller. 
             */
            unsigned int get_id() const;

        private: 
            state::StateManager* const owner;
            unsigned int id;
    };
}