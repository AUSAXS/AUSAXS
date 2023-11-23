#include <em/manager/ProteinManager.h>
#include <utility/Axis.h>
#include <em/ImageStack.h>
#include <settings/EMSettings.h>

using namespace em::managers;

ProteinManager::ProteinManager(std::observer_ptr<const em::ImageStackBase> images) : images(images) {
    double max = images->from_level(5);
    Axis axis(0, max, settings::em::charge_levels);
    set_charge_levels(axis.as_vector());
}

std::vector<double> ProteinManager::get_charge_levels() const noexcept {
    return charge_levels;
}

void ProteinManager::set_charge_levels(std::vector<double> levels) noexcept {
    // make sure the last bin can contain all atoms
    if (std::abs(levels[levels.size()-1]) < 10000) {
        levels.push_back(levels[0] < 0 ? -10000 : 10000);
    } 
    charge_levels = levels;
}