#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/RigidBodyManager.h>
#include <io/ExistingFile.h>

#include <iostream>

using namespace rigidbody::sequencer;

template<typename T> requires std::is_same_v<std::decay_t<T>, Protein>
Sequencer::Sequencer(const io::ExistingFile& saxs, T&& rigidbody) {
    rigidbody::sequencer::rigidbody = std::make_unique<RigidBodyManager>(saxs, std::forward<T>(rigidbody));
}
template Sequencer::Sequencer(const io::ExistingFile& saxs, Protein&& protein);
template Sequencer::Sequencer(const io::ExistingFile& saxs, const Protein& protein);
Sequencer::~Sequencer() = default;

void Sequencer::execute() {
    std::cout << "Sequencer::execute()" << std::endl;
    for (auto& loop : inner_loops) {
        loop->execute();
    }
}
