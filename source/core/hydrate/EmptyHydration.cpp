/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/EmptyHydration.h>

using namespace ausaxs;

hydrate::EmptyHydration::EmptyHydration() = default;

hydrate::EmptyHydration::~EmptyHydration() = default;

void hydrate::EmptyHydration::clear() {}

std::unique_ptr<hydrate::Hydration> hydrate::EmptyHydration::clone() const {return std::make_unique<EmptyHydration>();}