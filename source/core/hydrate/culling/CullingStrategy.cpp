/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/culling/CullingStrategy.h>

using namespace hydrate;

CullingStrategy::CullingStrategy(observer_ptr<data::Molecule> molecule) : molecule(molecule) {}

CullingStrategy::~CullingStrategy() = default;

void CullingStrategy::set_target_count(unsigned int target_count) {this->target_count = target_count;}
