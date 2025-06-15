#include <rigidbody/DefaultOptimizer.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/Rigidbody.h>
#include <fitter/FitResult.h>
#include <io/ExistingFile.h>
#include <settings/GeneralSettings.h>

using namespace ausaxs;

std::shared_ptr<fitter::FitResult> rigidbody::default_optimize(observer_ptr<Rigidbody> rigidbody, const io::ExistingFile& measurement_path) {
    rigidbody::sequencer::Sequencer sequencer(measurement_path);
    return sequencer
        .setup()
            .load_existing(rigidbody)
        .end()
        .save(settings::general::output + "initial.pdb")
        .loop(100)
            .optimize()
        .end()
        .save(settings::general::output + "optimized.pdb")
    .execute();
}