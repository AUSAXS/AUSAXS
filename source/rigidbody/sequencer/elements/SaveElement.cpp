// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/SaveElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/Rigidbody.h>
#include <settings/GeneralSettings.h>
#include <io/detail/trajectory/XYZWriter.h>

#include <unordered_map>

using namespace ausaxs::rigidbody::sequencer;

SaveElement::SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, const io::File& path) : LoopElementCallback(owner), path(path) {}
SaveElement::~SaveElement() = default;

void SaveElement::run() {
    if (const auto& ext = path.extension(); ext == ".pdb") {
        static int counter = 0;
        owner->_get_rigidbody()->molecule.save(path.append(std::to_string(counter++)));
    } else if (ext == ".xyz") {
        static std::unordered_map<std::string, io::detail::xyz::XYZWriter> writers;
        auto p = path.path(); 
        if (!writers.contains(p)) {
            writers.emplace(p, path);
        }
        writers.at(p).write_frame(&owner->_get_rigidbody()->molecule);
    } else {
        throw std::runtime_error("SaveElement::run: Unknown file format: \"" + ext + "\"");
    }
}
