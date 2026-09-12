// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <em/manager/ProteinManagerFactory.h>

#include <data/Molecule.h>  // IWYU pragma: keep
#include <em/manager/SmartProteinManager.h>

using namespace ausaxs;

std::unique_ptr<em::managers::ProteinManager> em::factory::create_manager(observer_ptr<const ImageStackBase> images) {
    // return std::make_unique<em::managers::SimpleProteinManager>(images);
    return std::make_unique<em::managers::SmartProteinManager>(images);
}