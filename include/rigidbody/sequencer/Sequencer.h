#pragma once

#include <rigidbody/sequencer/SequenceElement.h>

#include <memory>
#include <vector>

namespace rigidbody {
    namespace sequencer {
        class Sequencer {
            public:
                Sequencer() = default;
                virtual ~Sequencer() = default;

                Sequencer& add(std::shared_ptr<SequenceElement> element) {
                    return *this;
                }

                Sequencer& add(Loop loop) {
                    return this->loop(loop);
                }                

                Sequencer& loop(Loop loop) {
                    return *this;
                }                

                void execute();

            private: 
                std::vector<std::shared_ptr<SequenceElement>> elements;
        };
    }
}