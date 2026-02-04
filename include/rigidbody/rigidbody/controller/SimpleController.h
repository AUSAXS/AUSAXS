// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/controller/IController.h>

namespace ausaxs::rigidbody::controller {
    class SimpleController : public IController {
        public:
            using IController::IController;
            ~SimpleController() override;

            void setup(const io::ExistingFile& measurement_path) override;
            bool prepare_step() override;
            void finish_step() override;

        private:
            /**
             * @brief Update the fitter with the current histogram. 
             */
            void update_fitter(); 
    };
}