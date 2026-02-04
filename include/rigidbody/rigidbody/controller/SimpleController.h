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
            void prepare_step() override;
            bool finish_step() override;

        private:
            /**
             * @brief Update the fitter with the current histogram. 
             */
            void update_fitter(); 
    };
}