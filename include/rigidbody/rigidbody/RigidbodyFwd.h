#pragma once

namespace ausaxs::rigidbody {
    struct Rigidbody;
    namespace detail            {struct Configuration;}
    namespace transform         {class TransformStrategy;}
    namespace selection         {class BodySelectStrategy;}
    namespace constraints       {struct ConstraintManager;}
    namespace parameter         {class ParameterGenerationStrategy;}
    namespace parameter::decay  {class DecayStrategy;}
}

namespace ausaxs::fitter {
    struct ConstrainedFitter;
}