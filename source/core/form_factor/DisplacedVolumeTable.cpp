#include <form_factor/DisplacedVolumeTable.h>
#include <settings/MoleculeSettings.h>

#include <stdexcept>

using namespace ausaxs;

constants::displaced_volume::detail::DisplacedVolumeSet constants::displaced_volume::get_displaced_volume_set() {
    switch (settings::molecule::displaced_volume_set) {
        case settings::molecule::DisplacedVolumeSet::Traube: return Traube;
        case settings::molecule::DisplacedVolumeSet::Voronoi_explicit_H: return Voronoi_explicit_H;
        case settings::molecule::DisplacedVolumeSet::Voronoi_implicit_H: return Voronoi_implicit_H;
        case settings::molecule::DisplacedVolumeSet::MinimumFluctutation_explicit_H: return MinimumFluctuation_explicit_H;
        case settings::molecule::DisplacedVolumeSet::MinimumFluctutation_implicit_H: return MinimumFluctuation_implicit_H;
        case settings::molecule::DisplacedVolumeSet::vdw: return vdw;
        default: throw std::runtime_error("constants::displaced_volume::get_displaced_volume_set: Invalid displaced volume set (enum " + std::to_string(static_cast<int>(settings::molecule::displaced_volume_set)) + ")");
    }
}