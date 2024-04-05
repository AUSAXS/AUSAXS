/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/RigidBodyManager.h>
#include <io/ExistingFile.h>

#include <iostream>

using namespace rigidbody::sequencer;

template<typename T> requires std::is_base_of_v<data::Molecule, std::decay_t<T>>
Sequencer::Sequencer(const io::ExistingFile& saxs, T&& rigidbody) {
    rigidbody::sequencer::rigidbody = std::make_unique<RigidBodyManager>(saxs, std::forward<T>(rigidbody));
}
Sequencer::~Sequencer() = default;

std::shared_ptr<fitter::Fit> Sequencer::execute() {
    std::cout << "Sequencer::execute()" << std::endl;
    rigidbody->initialize();
    for (auto& loop : elements) {
        loop->run();
    }
    return rigidbody->get_fit();
}

template rigidbody::sequencer::Sequencer::Sequencer(const io::ExistingFile& saxs, data::Molecule&& protein);
template rigidbody::sequencer::Sequencer::Sequencer(const io::ExistingFile& saxs, const data::Molecule& protein);
template rigidbody::sequencer::Sequencer::Sequencer(const io::ExistingFile& saxs, data::Molecule& protein);
template rigidbody::sequencer::Sequencer::Sequencer(const io::ExistingFile& saxs, rigidbody::RigidBody&& protein);
template rigidbody::sequencer::Sequencer::Sequencer(const io::ExistingFile& saxs, const rigidbody::RigidBody& protein);
template rigidbody::sequencer::Sequencer::Sequencer(const io::ExistingFile& saxs, rigidbody::RigidBody& protein);