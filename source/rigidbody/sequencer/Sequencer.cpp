/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/RigidBodyManager.h>
#include <io/ExistingFile.h>

#include <iostream>

using namespace rigidbody::sequencer;

template<typename T> requires std::is_same_v<std::decay_t<T>, data::Molecule>
Sequencer::Sequencer(const io::ExistingFile& saxs, T&& rigidbody) {
    rigidbody::sequencer::rigidbody = std::make_unique<RigidBodyManager>(saxs, std::forward<T>(rigidbody));
}
template Sequencer::Sequencer(const io::ExistingFile& saxs, data::Molecule&& protein);
template Sequencer::Sequencer(const io::ExistingFile& saxs, const data::Molecule& protein);
Sequencer::~Sequencer() = default;

void Sequencer::execute() {
    std::cout << "Sequencer::execute()" << std::endl;
    for (auto& loop : inner_loops) {
        loop->execute();
    }
}
