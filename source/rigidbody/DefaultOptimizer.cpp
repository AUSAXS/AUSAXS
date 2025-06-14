#include <rigidbody/DefaultOptimizer.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <fitter/FitResult.h>
#include <io/ExistingFile.h>
#include <settings/GeneralSettings.h>

using namespace ausaxs;

std::shared_ptr<fitter::FitResult> rigidbody::default_optimize(observer_ptr<Rigidbody> rigidbody, const io::ExistingFile& measurement_path) {
    rigidbody->molecule.save(settings::general::output + "initial.pdb");
    auto sequencer = std::make_unique<rigidbody::sequencer::Sequencer>(measurement_path);
    sequencer->setup()->_set_active_body(rigidbody);
    sequencer->loop(1000).optimize();
    rigidbody->molecule.save(settings::general::output + "optimized.pdb");
    return sequencer->execute();
}