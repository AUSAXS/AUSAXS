/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/manager/ProteinManager.h>
#include <em/ImageStack.h>
#include <data/Molecule.h>
#include <utility/Axis.h>
#include <settings/EMSettings.h>
#include <settings/MoleculeSettings.h>

using namespace ausaxs;
using namespace ausaxs::em::managers;

ProteinManager::~ProteinManager() = default;

ProteinManager::ProteinManager(observer_ptr<const em::ImageStackBase> images) : images(images) {
    settings::molecule::center = false;                 // centering doesn't make sense for dummy structures
    settings::molecule::implicit_hydrogens = false;     // likewise we don't know how many hydrogens are attached
    double max = images->from_level(settings::em::alpha_levels.max);
    double min = images->from_level(settings::em::alpha_levels.min);
    Axis axis(min, max, settings::em::charge_levels);
    set_charge_levels(axis.as_vector());
}

double ProteinManager::get_volume_grid() const {
    auto protein = get_protein();
    assert(protein != nullptr && "ProteinManager::get_volume_grid: protein is null");
    return protein->get_volume_grid();
}

double ProteinManager::get_excluded_volume_mass() const {
    auto protein = get_protein();
    assert(protein != nullptr && "ProteinManager::get_excluded_volume_mass: protein is null");
    return protein->get_excluded_volume_mass();
}    

std::vector<double> ProteinManager::get_charge_levels() const noexcept {
    return charge_levels;
}

void ProteinManager::set_charge_levels(const std::vector<double>& levels) noexcept {
    auto tmp = levels;

    // make sure the last bin can contain all atoms
    if (std::abs(levels.back()) < 10000) {
        tmp.push_back(10000);
    } 
    charge_levels = std::move(tmp);
}