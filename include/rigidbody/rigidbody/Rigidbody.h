// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/Molecule.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/controller/IController.h>
#include <fitter/FitterFwd.h>

#include <memory>

namespace ausaxs::rigidbody {
	struct Rigidbody {
        Rigidbody(Rigidbody&& other);
        Rigidbody& operator=(Rigidbody&& other);
        Rigidbody(data::Molecule&& molecule);
        ~Rigidbody();

        data::Molecule molecule;
        std::unique_ptr<constraints::ConstraintManager> constraints;
        std::unique_ptr<controller::IController> controller;
        std::unique_ptr<detail::Conformation> conformation;

        // the following are shared because they may be owned by the sequencer
        std::shared_ptr<selection::BodySelectStrategy> body_selector;
        std::shared_ptr<transform::TransformStrategy> transformer;
        std::shared_ptr<parameter::ParameterGenerationStrategy> parameter_generator;

        /**
         * @brief Refresh the grid of the molecule if necessary. 
         *        This guarantees that there is room for the current conformation. 
         */
        void refresh_grid();
    };
}