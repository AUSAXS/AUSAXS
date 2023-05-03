#pragma once

#include <rigidbody/sequencer/SequenceElement.h>
#include <rigidbody/sequencer/LoopElement.h>

#include <memory>
#include <vector>

namespace rigidbody {
    namespace sequencer {
        class Sequencer : public LoopElement {
            public:
                Sequencer() : LoopElement(1) {
                    std::cout << "Sequencer::Sequencer()" << std::endl;
                };
                virtual ~Sequencer() = default;

                void execute() override {
                    std::cout << "Sequencer::execute()" << std::endl;
                    for (auto& loop : inner_loops) {
                        loop->execute();
                    }
                }
        };
    }
}