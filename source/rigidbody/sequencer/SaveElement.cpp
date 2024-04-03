#include <rigidbody/sequencer/SaveElement.h>
#include <settings/GeneralSettings.h>

using namespace rigidbody::sequencer;

SaveElement::SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner) : LoopElementCallback(owner) {};
SaveElement::SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, SaveFormat fmt) : LoopElementCallback(owner) {
    switch (fmt) {
        case SaveFormat::PDB: 
            write_pdb();
            break;
        case SaveFormat::XYZ:
            write_xyz();
            break;
    }
}
SaveElement::~SaveElement() override = default;

void write_pdb(const io::File& path);

void write_xyz(const io::File& path = "");