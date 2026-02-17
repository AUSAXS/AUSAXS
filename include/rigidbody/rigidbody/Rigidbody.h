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
        std::unique_ptr<detail::SystemSpecification> conformation;

        // the following are shared because they may be owned by the sequencer
        std::shared_ptr<selection::BodySelectStrategy> body_selector;
        std::shared_ptr<transform::TransformStrategy> transformer;
        std::shared_ptr<parameter::ParameterGenerationStrategy> parameter_generator;

        /**
         * @brief Refresh the grid of the molecule if necessary. 
         *        This guarantees that there is room for the current conformation. 
         */
        void refresh_grid();

        /**
         * @brief Generate the current state of a body from its centered initial conformation
         *        and absolute transformation parameters.
         * @param ibody Index of the body to generate.
         * @return Body in its current transformed state with symmetry applied.
         */
        data::Body generate_current_state(int ibody) const;

        /**
         * @brief Save the current optimized structure including symmetry transformations.
         *        This applies symmetry to the centered initial conformation and then applies
         *        the absolute transformation to get the actual optimized structure.
         */
        void save(const io::File& path) const;
    };
}