#pragma once

namespace rigidbody {
    class RigidBody;
    namespace detail{struct BestConf;}
    namespace transform {class TransformStrategy;}
    namespace selection {class BodySelectStrategy;}
    namespace constraints {class ConstraintManager;}
    namespace parameter {class ParameterGenerationStrategy;}
    namespace parameter::decay {class DecayStrategy;}
}