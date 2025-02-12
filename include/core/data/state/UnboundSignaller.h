#pragma once

#include <data/state/DataStateFwd.h>
#include <data/state/Signaller.h>

#include <iostream>
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

            void modified_internal() const override {}

            void modified_external() const override {}

            void modified_symmetry(int) const override {}

            void set_symmetry_size(std::size_t) const override {
                std::cout << "UnboundSignaller::set_symmetry_size: called on an unbound signaller." << std::endl;
            }
    };
}