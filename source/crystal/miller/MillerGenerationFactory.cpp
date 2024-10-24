/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <crystal/miller/MillerGenerationFactory.h>
#include <crystal/miller/AllMillers.h>
#include <crystal/miller/ReducedMillers.h>
#include <crystal/miller/FibonacciMillers.h>
#include <settings/CrystalSettings.h>
#include <residue/ResidueMap.h>

using namespace ausaxs;

std::unique_ptr<crystal::MillerGenerationStrategy> crystal::factory::construct_miller_strategy() {
    return construct_miller_strategy(settings::crystal::miller_generation_strategy);
}

std::unique_ptr<crystal::MillerGenerationStrategy> crystal::factory::construct_miller_strategy(const settings::crystal::MillerGenerationChoice& choice) {
    switch (choice) {
        case settings::crystal::MillerGenerationChoice::All: {
            return std::make_unique<crystal::AllMillers>(settings::crystal::h, settings::crystal::k, settings::crystal::l);
        }
        case settings::crystal::MillerGenerationChoice::Fibonacci: {
            return std::make_unique<crystal::FibonacciMillers>(settings::crystal::h, settings::crystal::k, settings::crystal::l);
        }
        case settings::crystal::MillerGenerationChoice::Reduced: {
            return std::make_unique<crystal::ReducedMillers>(settings::crystal::h, settings::crystal::k, settings::crystal::l);
        }
        default: {
            throw std::invalid_argument("crystal::factory::construct_miller_strategy: Unknown MillerGenerationStrategy. Did you forget to add it to the switch statement?");
        }
    }
}