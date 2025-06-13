#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>
#include <fitter/FitResult.h>
#include <io/IOFwd.h>

#include <memory>

namespace ausaxs::rigidbody {
    std::shared_ptr<fitter::FitResult> default_optimize(observer_ptr<RigidBody> rigidbody, const io::ExistingFile& measurement_path);
}