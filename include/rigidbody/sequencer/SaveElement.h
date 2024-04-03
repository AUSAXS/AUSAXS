#pragma once

#include <io/File.h>
#include <utility/observer_ptr.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>

namespace rigidbody {
    namespace sequencer {
        enum class SaveFormat {XYZ, PDB};

        class SaveElement : public LoopElementCallback {
            public:
                SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner);
                SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, SaveFormat fmt);
                ~SaveElement() override;

                /**
                 * @brief Save the current structure as a pdb file. 
                 * 
                 * @arg path The save location. For ease of use a unique number will automatically be added as a suffix. 
                 *           Leave empty to save at a default location. 
                 */
                void write_pdb(const io::File& path = "");

                /**
                 * @brief Append the current structure to the animated xyz file. 
                 * 
                 * @arg path The save location. If the file already exists, the current state will be added. 
                 *           Leave empty to save at a default location. 
                 */
                void write_xyz(const io::File& path = "");
        };
    }
}