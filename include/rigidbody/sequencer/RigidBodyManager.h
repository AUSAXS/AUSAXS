#pragma once

#include <io/IOFwd.h>
#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <rigidbody/RigidBody.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace rigidbody::sequencer {
    class RigidBodyManager : public RigidBody {
        public:
            template<typename T> requires std::is_base_of_v<data::Molecule, std::decay_t<T>>
            RigidBodyManager(const io::ExistingFile& saxs, T&& rigidbody);
            ~RigidBodyManager();

            void optimize_step();

            void set_constraint_manager(std::shared_ptr<rigidbody::constraints::ConstraintManager> constraints);

            void set_body_select_manager(std::shared_ptr<rigidbody::selection::BodySelectStrategy> body_selector);

            void set_transform_manager(std::shared_ptr<rigidbody::transform::TransformStrategy> transform);

            void set_parameter_manager(std::shared_ptr<rigidbody::parameter::ParameterGenerationStrategy> parameters);

            void initialize();

            std::shared_ptr<fitter::Fit> get_fit() const;

        private:
            detail::BestConf best;
    };
    extern std::unique_ptr<RigidBodyManager> rigidbody;
}