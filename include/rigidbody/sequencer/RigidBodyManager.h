#pragma once

#include <rigidbody/RigidBody.h>
#include <io/ExistingFile.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace rigidbody {
    class RigidBody;
    namespace sequencer {
        class RigidBodyManager : RigidBody {
            public:
                template<typename T> requires std::is_same_v<std::decay_t<T>, Protein>
                RigidBodyManager(const io::ExistingFile& saxs, T&& rigidbody);
                ~RigidBodyManager();

                void optimize_step();

                void set_managers(settings::rigidbody::BodySelectStrategyChoice body_selector, 
                                    settings::rigidbody::TransformationStrategyChoice transform,
                                    settings::rigidbody::ParameterGenerationStrategyChoice parameters
                );

            private:
                rigidbody::detail::BestConf best;
                void initialize();
        };
        std::unique_ptr<RigidBodyManager> rigidbody;
    }
}