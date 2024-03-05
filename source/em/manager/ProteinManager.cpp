#include <em/manager/ProteinManager.h>
#include <utility/Axis.h>
#include <em/ImageStack.h>
#include <settings/EMSettings.h>
#include <settings/MoleculeSettings.h>

using namespace em::managers;

ProteinManager::~ProteinManager() = default;

ProteinManager::ProteinManager(observer_ptr<const em::ImageStackBase> images) : images(images) {
    settings::molecule::center = false;                 // centering doesn't make sense for dummy structures
    settings::molecule::use_effective_charge = false;   // we don't know the actual charge of a dummy structure
    settings::molecule::implicit_hydrogens = false;     // likewise we don't know how many hydrogens are attached
    double max = images->from_level(settings::em::alpha_levels.max);
    double min = images->from_level(settings::em::alpha_levels.min);
    Axis axis(min, max, settings::em::charge_levels);
    set_charge_levels(axis.as_vector());
}

std::vector<double> ProteinManager::get_charge_levels() const noexcept {
    return charge_levels;
}

void ProteinManager::set_charge_levels(const std::vector<double>& levels) noexcept {
    auto tmp = levels;

    // make sure the last bin can contain all atoms
    if (std::abs(levels.back()) < 10000) {
        tmp.push_back(levels.front() < 0 ? -10000 : 10000);
    } 
    charge_levels = tmp;
}