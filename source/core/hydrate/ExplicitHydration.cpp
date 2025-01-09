/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include "hydrate/Hydration.h"
#include <hydrate/ExplicitHydration.h>
#include <data/atoms/Water.h>

using namespace ausaxs;

hydrate::ExplicitHydration::ExplicitHydration() = default;
hydrate::ExplicitHydration::ExplicitHydration(const std::vector<data::Water>& waters) : waters(waters) {}
hydrate::ExplicitHydration::ExplicitHydration(std::vector<data::Water>&& waters) : waters(std::move(waters)) {}
hydrate::ExplicitHydration::~ExplicitHydration() = default;

void hydrate::ExplicitHydration::clear() {
    waters.clear();
}

std::unique_ptr<hydrate::Hydration> hydrate::ExplicitHydration::clone() const {
    return std::make_unique<hydrate::ExplicitHydration>(*this);
}