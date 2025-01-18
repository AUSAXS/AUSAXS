#pragma once

#include <data/state/DataStateFwd.h>
#include <data/state/Signaller.h>

namespace ausaxs::signaller {
    /**
     * @brief Dummy version of a Signaller object. This can be used to initialize an instance of Signaller. 
     */
    class UnboundSignaller : public Signaller {
        public: 
			UnboundSignaller() = default;
			UnboundSignaller(const UnboundSignaller& rhs) = default;
			UnboundSignaller(UnboundSignaller&& rhs) noexcept = default;
			UnboundSignaller &operator=(const UnboundSignaller& rhs) = default;
			UnboundSignaller &operator=(UnboundSignaller&& rhs) noexcept = default;
			~UnboundSignaller() override = default;

            void internal_change() const override {}

            void external_change() const override {}

            void symmetry_changed() const override {}
    };
}