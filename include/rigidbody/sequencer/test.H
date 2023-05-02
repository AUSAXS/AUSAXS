#pragma once

#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/Element.h>
#include <rigidbody/sequencer/Sequencer.h>

void test() {
    rigidbody::sequencer::Sequencer()
        .loop(rigidbody::sequencer::Loop()
            .add(std::make_shared<rigidbody::sequencer::Parameters>())
            .repeat(10)
        )
    .execute();
}