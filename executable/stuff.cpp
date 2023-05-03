#include <rigidbody/sequencer/Sequencer.h>

int main(int argc, char const *argv[]) {
    rigidbody::sequencer::Sequencer()
        .loop(5)
            .body_select_strategy(settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect)
            .parameter_strategy(settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly)
            .transform_strategy(settings::rigidbody::TransformationStrategyChoice::RigidTransform)
    .execute();

    return 0;
}