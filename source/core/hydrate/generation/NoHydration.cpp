/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/generation/NoHydration.h>
#include <hydrate/ExplicitHydration.h>
#include <data/record/Water.h>

using namespace ausaxs;

std::unique_ptr<hydrate::Hydration> hydrate::NoHydration::hydrate() {
    return std::make_unique<ExplicitHydration>(std::vector<data::record::Water>());
}