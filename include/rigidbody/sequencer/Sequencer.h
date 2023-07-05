#pragma once

#include <rigidbody/sequencer/LoopElement.h>
#include <data/Protein.h>

#include <memory>
#include <vector>
#include <concepts>

namespace rigidbody {
    namespace sequencer {
        class Sequencer : public LoopElement {
            public:
                template<typename T> requires std::is_same_v<std::decay_t<T>, Protein>
                Sequencer(const io::ExistingFile& saxs, T&& Protein);
                ~Sequencer();

                void execute() override;
        };
    }
}