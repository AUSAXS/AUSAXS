#pragma once

namespace ausaxs::rigidbody {
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
        struct SymmetryParameter;
    }

    namespace constraints {
        class ConstraintManager;
        class DistanceConstraint;
    }
}