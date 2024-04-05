#pragma once

#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/sequencer/ConstraintIteratorElement.h>
#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/TransformElement.h>
#include <data/Molecule.h>

namespace rigidbody {
    namespace sequencer {
        class Sequencer : public LoopElement {
            public:
                template<typename T> requires std::is_base_of_v<data::Molecule, std::decay_t<T>>
                Sequencer(const io::ExistingFile& saxs, T&& protein);
                ~Sequencer();

                std::shared_ptr<fitter::Fit> execute() override;
        };
    }
}