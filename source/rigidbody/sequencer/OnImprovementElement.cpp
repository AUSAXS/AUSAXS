/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/OnImprovementElement.h>
#include <rigidbody/detail/BestConf.h>

#include <limits>

using namespace ausaxs::rigidbody::sequencer;

OnImprovementElement::OnImprovementElement(observer_ptr<LoopElement> owner) : LoopElement(owner, 1), best_chi2(std::numeric_limits<double>::max()) {}

OnImprovementElement::~OnImprovementElement() = default;

void OnImprovementElement::run() {
    double new_best_conf = _get_best_conf()->chi2;
    if (new_best_conf < best_chi2) {
        best_chi2 = new_best_conf;
        for (auto& e : elements) {
            e->run();
        }
    }
}