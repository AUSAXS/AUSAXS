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
SaveElement::~SaveElement() = default;

void SaveElement::write_pdb(const io::File& path) {
    if (path.empty()) {
        // save at default location
    } else {
        // save at specified location
    }
}

void SaveElement::write_xyz(const io::File& path = "") {

}