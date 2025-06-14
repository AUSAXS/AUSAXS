// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::rigidbody {
    struct Rigidbody;
    namespace detail            {struct BestConf;}
    namespace transform         {class TransformStrategy;}
    namespace selection         {class BodySelectStrategy;}
    namespace constraints       {class ConstraintManager;}
    namespace parameter         {class ParameterGenerationStrategy;}
    namespace parameter::decay  {class DecayStrategy;}
}

namespace ausaxs::fitter {
    class ConstrainedFitter;
}