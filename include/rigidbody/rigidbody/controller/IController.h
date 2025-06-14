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
            * @brief Run a single optimization step. 
            */
            virtual bool run_step() = 0;

            observer_ptr<detail::BestConf> current_best() const;
            observer_ptr<fitter::ConstrainedFitter> get_fitter() const;
            observer_ptr<const fitter::FitResult> get_calibration() const;

            protected:
                observer_ptr<Rigidbody> rigidbody;
                std::unique_ptr<fitter::ConstrainedFitter> fitter;
                std::unique_ptr<fitter::FitResult> calibration;
                std::unique_ptr<rigidbody::detail::BestConf> best;
    };
}