// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <fitter/FitterFwd.h>
#include <io/IOFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::rigidbody::controller {
    /**
     * @brief Interface for a controller that can be used to control the rigid body optimization process.
     */
    class IController {
        public:
            IController(observer_ptr<Rigidbody> rigidbody);
            virtual ~IController();

            /**
            * @brief Setup the controller.
            * 
            * This method will be called before the optimization process starts, allowing the controller to perform any necessary setup.
            */
            virtual void setup(const io::ExistingFile& measurement_path) = 0;

            /**
            * @brief Prepare the next optimization step. This prepares the internal state for the next optimization step.
            *
            * The step is split in two to provide access to the generated configuration before it is evaluated.
            * finish_step() must be called to complete the step.
            * @return true if the step will be accepted, false otherwise.
            */
            virtual bool prepare_step() = 0;

            /**
             * @brief Finish the current step. 
             */
            virtual void finish_step() = 0;

            /**
             * @brief Load the histogram of the molecule's current state into the fitter.
             *
             * The fitter otherwise holds the model of the last evaluated candidate, which after a rejected step is not
             * the state the molecule was restored to. Cheap when the molecule has not changed since the last update.
             */
            virtual void update_fitter() = 0;

            /// @brief Get the best body configuration found so far.
            observer_ptr<detail::MoleculeTransformParametersAbsolute> get_current_best_config() const;

            /// @brief Get the fitter used to evaluate candidate configurations.
            observer_ptr<fitter::ConstrainedFitter> get_fitter() const;

        protected:
            bool step_accepted = false;
            observer_ptr<Rigidbody> rigidbody;
            std::unique_ptr<fitter::ConstrainedFitter> fitter;
            std::unique_ptr<rigidbody::detail::MoleculeTransformParametersAbsolute> current_best_config;            
    };
}