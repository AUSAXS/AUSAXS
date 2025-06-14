#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>
#include <fitter/FitterFwd.h>
#include <io/IOFwd.h>

#include <memory>

namespace ausaxs::rigidbody {
    std::shared_ptr<fitter::FitResult> default_optimize(observer_ptr<Rigidbody> rigidbody, const io::ExistingFile& measurement_path);
}