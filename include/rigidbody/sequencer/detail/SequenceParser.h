#pragma once

#include <rigidbody/sequencer/Sequencer.h>
#include <io/IOFwd.h>

namespace rigidbody::sequencer {
    class SequenceParser {
        public:
            SequenceParser() = default;

            Sequencer parse(const io::ExistingFile& file);
    };
}