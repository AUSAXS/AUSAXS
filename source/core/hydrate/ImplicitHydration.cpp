/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/ImplicitHydration.h>

using namespace ausaxs;

hydrate::ImplicitHydration::ImplicitHydration() = default;

hydrate::ImplicitHydration::~ImplicitHydration() = default;

void hydrate::ImplicitHydration::clear() {throw std::runtime_error("ImplicitHydration::clear: Not implemented.");}

std::unique_ptr<hydrate::Hydration> hydrate::ImplicitHydration::clone() const {throw std::runtime_error("ImplicitHydration::clone: Not implemented.");}