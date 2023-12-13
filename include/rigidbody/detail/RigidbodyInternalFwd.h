#pragma once

namespace rigidbody {
    namespace detail {
        struct BestConf;
    }

	namespace selection {
        class BodySelectStrategy;
    }

    namespace transform {
    	class TransformStrategy;
        struct TransformGroup;
        struct BackupBody;
    }

    namespace parameter {
    	class ParameterGenerationStrategy;
    }

    namespace constraints {
        class ConstraintManager;
        class DistanceConstraint;
    }
}