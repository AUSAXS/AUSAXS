/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/SaveElement.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/RigidBody.h>
#include <settings/GeneralSettings.h>
#include <io/XYZWriter.h>

#include <unordered_map>

using namespace rigidbody::sequencer;

SaveElement::SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, const io::File& path) : LoopElementCallback(owner), path(path) {}
SaveElement::~SaveElement() = default;

void SaveElement::run() {
    static int counter = 0;
    static std::unordered_map<std::string, io::XYZWriter> writers;
    if (const auto& ext = path.extension(); ext == ".pdb") {
        owner->_get_rigidbody()->save(path.append(std::to_string(counter++)));
    } else if (ext == ".xyz") {
        auto p = path.path(); 
        if (!writers.contains(p)) {
            writers.emplace(p, path);
        }
        writers.at(p).write_frame(owner->_get_rigidbody());
    } else {
        throw std::runtime_error("SaveElement::run: Unknown file format: \"" + ext + "\"");
    }
}
