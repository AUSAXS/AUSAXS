#include <rigidbody/sequencer/OptimizeStepElement.h>
#include <rigidbody/sequencer/RigidBodyManager.h>

#include <iostream>

namespace rigidbody::sequencer {
    OptimizeStepElement::OptimizeStepElement(LoopElement* owner) : LoopElementCallback(owner) {}

    OptimizeStepElement::~OptimizeStepElement() = default;

    void OptimizeStepElement::run() {
        std::cout << "OptimizeStepElement::run()" << std::endl;
        rigidbody->optimize_step();
    }
}