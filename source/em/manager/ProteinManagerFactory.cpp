/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/manager/ProteinManagerFactory.h>
#include <em/manager/SmartProteinManager.h>
#include <em/manager/SimpleProteinManager.h>
#include <data/Molecule.h>

using namespace ausaxs;

std::unique_ptr<em::managers::ProteinManager> em::factory::create_manager(observer_ptr<const ImageStackBase> images) {
    return std::make_unique<em::managers::SimpleProteinManager>(images);
    // return std::make_unique<em::managers::SmartProteinManager>(images);
}