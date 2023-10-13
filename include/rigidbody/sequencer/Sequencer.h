#pragma once

#include <rigidbody/sequencer/LoopElement.h>
#include <data/Molecule.h>

#include <memory>
#include <vector>
#include <concepts>

namespace rigidbody {
    namespace sequencer {
        class Sequencer : public LoopElement {
            public:
                template<typename T> requires std::is_same_v<std::decay_t<T>, data::Molecule>
                Sequencer(const io::ExistingFile& saxs, T&& protein);
                ~Sequencer();

                void execute() override;
        };
    }
}