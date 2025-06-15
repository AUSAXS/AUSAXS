#include <rigidbody/controller/ControllerFactory.h>
#include <rigidbody/controller/SimpleController.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;

std::unique_ptr<controller::IController> factory::create_controller(observer_ptr<Rigidbody> molecule) {
    return create_controller(molecule, settings::rigidbody::controller_choice);
}

std::unique_ptr<controller::IController> factory::create_controller(observer_ptr<Rigidbody> molecule, settings::rigidbody::ControllerChoice choice) {
    switch (choice) {
        case settings::rigidbody::ControllerChoice::Classic:
            return std::make_unique<controller::SimpleController>(molecule);
        case settings::rigidbody::ControllerChoice::Metropolis:
            throw std::runtime_error("rigidbody::factory::create_controller: Metropolis controller not implemented yet.");
        default: 
            throw except::unknown_argument("rigidbody::factory::create_controller: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}