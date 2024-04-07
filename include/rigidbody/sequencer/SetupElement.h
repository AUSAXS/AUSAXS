#pragma once

#include <rigidbody/sequencer/GenericElement.h>

namespace rigidbody::sequencer {
    class SetupElement : public GenericElement {
        public:
            SetupElement() = default;
            virtual ~SetupElement() = default;

            void run() override;
    };
}