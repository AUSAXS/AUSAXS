/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/state/UnboundSignaller.h>

using namespace signaller;

UnboundSignaller::UnboundSignaller() = default;

UnboundSignaller::~UnboundSignaller() = default;

void UnboundSignaller::external_change() const {}

void UnboundSignaller::internal_change() const {}
