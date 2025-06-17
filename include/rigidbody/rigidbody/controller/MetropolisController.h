#pragma once

#include <rigidbody/controller/IController.h>

namespace ausaxs::rigidbody::controller {
    class MetropolisController : public IController {
        public:
            using IController::IController;
            ~MetropolisController() override;

            void setup(const io::ExistingFile& measurement_path) override; //< @copydoc IController::setup()
            bool run_step() override; //< @copydoc IController::run_step()

        private:
            /**
             * @brief Update the fitter with the current histogram. 
             */
            void update_fitter(); 
    };
}