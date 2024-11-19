#include <hist/detail/SimpleExvModel.h>
#include <hist/detail/CompactCoordinates.h>
#include <data/Molecule.h>

using namespace ausaxs::hist::detail;

SimpleExvModel::SimpleExvModel(observer_ptr<const data::Molecule> protein) : protein(protein) {}
SimpleExvModel::~SimpleExvModel() = default;

bool flag_simple_excluded_volume = false;
void SimpleExvModel::enable() {
    flag_simple_excluded_volume = true;
}

void SimpleExvModel::disable() {
    flag_simple_excluded_volume = false;
}

void SimpleExvModel::apply_simple_excluded_volume(hist::detail::CompactCoordinates& data_a) {
    if (flag_simple_excluded_volume) {
        data_a.implicit_excluded_volume(protein->get_volume_grid()/protein->size_atom());
    }
}