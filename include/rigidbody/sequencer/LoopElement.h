#pragma once

#include <rigidbody/sequencer/SequenceElement.h>

#include <memory>
#include <vector>

namespace rigidbody {
    namespace sequencer {
        class Loop {
            public:
                Loop() = default;
                virtual ~Loop() = default;

                void execute() {
                    for (auto& element : elements) {
                        element->execute();
                    }
                }

                Loop& add(std::shared_ptr<SequenceElement> element) {
                    elements.push_back(element);
                    return *this;
                }

                Loop& repeat(unsigned int iterations) {
                    this->iterations = iterations;
                    return *this;
                }

            private: 
                std::vector<std::shared_ptr<SequenceElement>> elements;
                unsigned int iterations = 1;
        };
    }
}