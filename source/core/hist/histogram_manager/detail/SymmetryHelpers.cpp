// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/detail/SymmetryHelpers.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::symmetry::detail;
using namespace ausaxs::hist::detail;

template<bool variable_bin_width>
std::pair<std::vector<BodySymmetryData<variable_bin_width>>, hist::detail::CompactCoordinates<variable_bin_width>> ausaxs::symmetry::detail::generate_transformed_data(const data::Molecule& protein) {
    std::vector<BodySymmetryData<variable_bin_width>> res(protein.size_body());
    for (int i_body1 = 0; i_body1 < static_cast<int>(protein.size_body()); ++i_body1) {
        res[i_body1] = generate_transformed_data<variable_bin_width>(protein.get_body(i_body1));
    }
    return {std::move(res), protein.get_waters()};
}
template std::pair<std::vector<BodySymmetryData<true>>, hist::detail::CompactCoordinates<true>> ausaxs::symmetry::detail::generate_transformed_data<true>(const data::Molecule& protein);
template std::pair<std::vector<BodySymmetryData<false>>, hist::detail::CompactCoordinates<false>> ausaxs::symmetry::detail::generate_transformed_data<false>(const data::Molecule& protein);

template<bool variable_bin_width>
BodySymmetryData<variable_bin_width> ausaxs::symmetry::detail::generate_transformed_data(const data::Body& body) {
    CompactCoordinates<variable_bin_width> data_a(body.get_atoms());
    auto cm = body.get_cm();

    // loop over its symmetries
    std::vector<std::vector<CompactCoordinates<variable_bin_width>>> atomic(1+body.size_symmetry());
    for (int i_sym_1 = 0; i_sym_1 < static_cast<int>(body.size_symmetry()); ++i_sym_1) {
        auto symmetry = body.symmetry().get(i_sym_1);

        // for every symmetry, loop over how many times it should be repeated
        // it is then repeatedly applied to the same data
        std::vector<CompactCoordinates<variable_bin_width>> sym_atomic(static_cast<int>(symmetry->repetitions()), data_a);
        for (int i_repeat = 0; i_repeat < static_cast<int>(symmetry->repetitions()); ++i_repeat) {
            auto t = symmetry->get_transform(cm, i_repeat+1);
            std::transform(
                sym_atomic[i_repeat].get_data().begin(), 
                sym_atomic[i_repeat].get_data().end(), 
                sym_atomic[i_repeat].get_data().begin(), 
                [t=std::move(t)] (const CompactCoordinatesXYZW<variable_bin_width>& v) -> CompactCoordinatesXYZW<variable_bin_width> {return {t(v.value.pos), v.value.w}; }
            );
        }
        atomic[1+i_sym_1] = std::move(sym_atomic);
    }

    atomic[0] = {std::move(data_a)};
    return {std::move(atomic)};
}
template BodySymmetryData<true> ausaxs::symmetry::detail::generate_transformed_data<true>(const data::Body& body);
template BodySymmetryData<false> ausaxs::symmetry::detail::generate_transformed_data<false>(const data::Body& body);

template<bool variable_bin_width>
SymmetryData<variable_bin_width> ausaxs::symmetry::detail::generate_transformed_data(const data::Body& body, int isym) {
    CompactCoordinates<variable_bin_width> data_a(body.get_atoms());
    auto cm = body.get_cm();
    auto symmetry = body.symmetry().get(isym);

    std::vector<CompactCoordinates<variable_bin_width>> sym_atomic(static_cast<int>(symmetry->repetitions()), data_a);
    for (int i_repeat = 0; i_repeat < static_cast<int>(symmetry->repetitions()); ++i_repeat) {
        auto t = symmetry->get_transform(cm, i_repeat+1);
        std::transform(
            sym_atomic[i_repeat].get_data().begin(), 
            sym_atomic[i_repeat].get_data().end(), 
            sym_atomic[i_repeat].get_data().begin(), 
            [t=std::move(t)] (const CompactCoordinatesXYZW<variable_bin_width>& v) -> CompactCoordinatesXYZW<variable_bin_width> {return {t(v.value.pos), v.value.w}; }
        );
    }
    return {std::move(sym_atomic)}; 
}
template SymmetryData<true> ausaxs::symmetry::detail::generate_transformed_data<true>(const data::Body& body, int isym);
template SymmetryData<false> ausaxs::symmetry::detail::generate_transformed_data<false>(const data::Body& body, int isym);

template struct ausaxs::symmetry::detail::BodySymmetryData<true>;
template struct ausaxs::symmetry::detail::BodySymmetryData<false>;
template struct ausaxs::symmetry::detail::SymmetryData<true>;
template struct ausaxs::symmetry::detail::SymmetryData<false>;
