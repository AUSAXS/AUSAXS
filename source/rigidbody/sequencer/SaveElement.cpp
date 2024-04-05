#include <rigidbody/sequencer/SaveElement.h>
#include <rigidbody/sequencer/RigidBodyManager.h>
#include <settings/GeneralSettings.h>
#include <io/XYZWriter.h>

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
    static int counter = 0;
    if (path.path().empty()) {
        io::File path = settings::general::output + "models/" + std::to_string(counter) + ".pdb";
        rigidbody->save(path);
    } else {
        rigidbody->save(path);
    }
}

void SaveElement::write_xyz() {
    static io::XYZWriter writer(settings::general::output + "animated.xyz");
    writer.write_frame(rigidbody.get());
}