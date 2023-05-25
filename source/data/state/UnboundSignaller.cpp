#include <data/state/UnboundSignaller.h>

using namespace signaller;

UnboundSignaller::UnboundSignaller() = default;

UnboundSignaller::~UnboundSignaller() = default;

#include <iostream>
void UnboundSignaller::external_change() const {
    std::cout << "UNBOUNDSIGNALLER: External change" << std::endl;
}

void UnboundSignaller::internal_change() const {
    std::cout << "UNBOUNDSIGNALLER: Internal change" << std::endl;
}
