#pragma once

#include <data/Molecule.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/controller/IController.h>
#include <fitter/FitterFwd.h>

#include <memory>

namespace ausaxs::rigidbody {
	struct Rigidbody {
        Rigidbody(data::Molecule&& molecule);
        ~Rigidbody();

        data::Molecule molecule;
        std::unique_ptr<constraints::ConstraintManager> constraints;
        std::unique_ptr<controller::IController> controller;

        // the following are shared because they may be owned by the sequencer
        std::shared_ptr<selection::BodySelectStrategy> body_selector;
        std::shared_ptr<transform::TransformStrategy> transformer;
        std::shared_ptr<parameter::ParameterGenerationStrategy> parameter_generator;
    };
}