// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

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
        struct BodyTransformParameters;
    }

    namespace constraints {
        struct ConstraintManager;
        class DistanceConstraint;
    }
}