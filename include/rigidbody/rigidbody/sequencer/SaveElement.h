#pragma once

#include <io/File.h>
#include <utility/observer_ptr.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>

namespace rigidbody::sequencer {
    class SaveElement : public LoopElementCallback, public GenericElement {
        public:
            SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, const io::File& path);
            ~SaveElement() override;

            void run() override;

        private:
            io::File path;
    };
}