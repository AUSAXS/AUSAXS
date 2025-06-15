#pragma once

#include <rigidbody/controller/IController.h>
#include <fitter/FitterFwd.h>

namespace ausaxs::rigidbody::controller {
    class SimpleController : public IController {
        public:
            using IController::IController;
            ~SimpleController() override;

            void setup(const io::ExistingFile& measurement_path) override; //< @copydoc IController::setup()
            bool run_step() override; //< @copydoc IController::run_step()

        private:
            void update_fitter(); 
    };
}