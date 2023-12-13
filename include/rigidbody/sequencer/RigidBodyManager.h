#pragma once

#include <io/IOFwd.h>
#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <rigidbody/RigidBody.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace rigidbody {
    namespace sequencer {
        class RigidBodyManager : RigidBody {
            public:
                template<typename T> requires std::is_same_v<std::decay_t<T>, data::Molecule>
                RigidBodyManager(const io::ExistingFile& saxs, T&& rigidbody);
                ~RigidBodyManager();

                void optimize_step();

                void set_managers(const settings::rigidbody::BodySelectStrategyChoice& body_selector, 
                                    const settings::rigidbody::TransformationStrategyChoice& transform,
                                    const settings::rigidbody::ParameterGenerationStrategyChoice& parameters
                );

            private:
                detail::BestConf best;
                void initialize();
        };
        extern std::unique_ptr<RigidBodyManager> rigidbody;
    }
}