// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>
#include <fitter/FitterFwd.h>
#include <io/IOFwd.h>

#include <memory>

namespace ausaxs::rigidbody::controller {
    /**
     * @brief Interface for a controller that can be used to control the rigid body optimization process.
     */
    class IController {
        public:
            IController(observer_ptr<Rigidbody> rigidbody);
            IController(observer_ptr<Rigidbody> rigidbody, std::unique_ptr<fitter::FitResult> calibration);
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

            observer_ptr<detail::Configuration> get_current_best_config() const;
            observer_ptr<detail::Configuration> get_current_config() const;
            observer_ptr<fitter::ConstrainedFitter> get_fitter() const;
            observer_ptr<const fitter::FitResult> get_calibration() const;

            protected:
                bool step_accepted = false;
                observer_ptr<Rigidbody> rigidbody;
                std::unique_ptr<fitter::ConstrainedFitter> fitter;
                std::unique_ptr<fitter::FitResult> calibration;
                std::unique_ptr<rigidbody::detail::Configuration> current_config;
                std::unique_ptr<rigidbody::detail::Configuration> current_best_config;            
    };
}