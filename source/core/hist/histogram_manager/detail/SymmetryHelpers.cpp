// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/detail/SymmetryHelpers.h>

#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/CompactCoordinatesFactoryFF.h>

#include <utility>

using namespace ausaxs;
using namespace ausaxs::symmetry::detail;
using namespace ausaxs::hist::detail;

namespace {
    template<bool form_factors>
    AtomicCoordinates<form_factors> construct(const data::Body& body) {
        if constexpr (form_factors) {return factory::construct_by_ff(body.get_atoms());}
        else                        {return factory::construct(body.get_atoms());}
    }

    void transform_all(CompactCoordinates& atoms, const transform::Affine& t) {atoms.transform_coordinates(t);}
    void transform_all(std::vector<CompactCoordinates>& atoms, const transform::Affine& t) {
        for (auto& set : atoms) {set.transform_coordinates(t);}
    }
}

template<bool form_factors>
std::pair<std::vector<BodySymmetryData<form_factors>>, hist::detail::CompactCoordinates> ausaxs::symmetry::detail::generate_transformed_data(const data::Molecule& protein) {
    std::vector<BodySymmetryData<form_factors>> res(protein.size_body());
    for (int i_body1 = 0; i_body1 < protein.size_body(); ++i_body1) {
        res[i_body1] = generate_transformed_data<form_factors>(protein.get_body(i_body1));
    }
    return {std::move(res), hist::detail::factory::construct_from_waters(&protein)};
}

template<bool form_factors>
BodySymmetryData<form_factors> ausaxs::symmetry::detail::generate_transformed_data(const data::Body& body) {
    auto data_a = construct<form_factors>(body);
    auto cm = body.get_cm();

    // loop over its symmetries
    std::vector<std::vector<AtomicCoordinates<form_factors>>> atomic(1+body.size_symmetry());
    for (int i_sym_1 = 0; i_sym_1 < body.size_symmetry(); ++i_sym_1) {
        const auto* symmetry = body.symmetry().get(i_sym_1);

        // for every symmetry, loop over how many times it should be repeated
        // it is then repeatedly applied to the same data
        std::vector<AtomicCoordinates<form_factors>> sym_atomic(symmetry->repetitions(), data_a);
        for (int i_repeat = 0; i_repeat < symmetry->repetitions(); ++i_repeat) {
            auto t = body.symmetry().get_transform(i_sym_1, cm, i_repeat+1);
            transform_all(sym_atomic[i_repeat], t);
        }
        atomic[1+i_sym_1] = std::move(sym_atomic);
    }

    atomic[0] = {std::move(data_a)};
    return {std::move(atomic)};
}

template<bool form_factors>
SymmetryData<form_factors> ausaxs::symmetry::detail::generate_transformed_data(const data::Body& body, int isym) {
    auto data_a = construct<form_factors>(body);
    auto cm = body.get_cm();
    const auto* symmetry = body.symmetry().get(isym);

    std::vector<AtomicCoordinates<form_factors>> sym_atomic(symmetry->repetitions(), data_a);
    for (int i_repeat = 0; i_repeat < symmetry->repetitions(); ++i_repeat) {
        auto t = body.symmetry().get_transform(isym, cm, i_repeat+1);
        transform_all(sym_atomic[i_repeat], t);
    }
    return {std::move(sym_atomic)};
}

template std::pair<std::vector<BodySymmetryData<false>>, hist::detail::CompactCoordinates> ausaxs::symmetry::detail::generate_transformed_data<false>(const data::Molecule&);
template BodySymmetryData<false> ausaxs::symmetry::detail::generate_transformed_data<false>(const data::Body&);
template SymmetryData<false> ausaxs::symmetry::detail::generate_transformed_data<false>(const data::Body&, int);
template struct ausaxs::symmetry::detail::BodySymmetryData<false>;
template struct ausaxs::symmetry::detail::SymmetryData<false>;
template std::pair<std::vector<BodySymmetryData<true>>, hist::detail::CompactCoordinates> ausaxs::symmetry::detail::generate_transformed_data<true>(const data::Molecule&);
template BodySymmetryData<true> ausaxs::symmetry::detail::generate_transformed_data<true>(const data::Body&);
template SymmetryData<true> ausaxs::symmetry::detail::generate_transformed_data<true>(const data::Body&, int);
template struct ausaxs::symmetry::detail::BodySymmetryData<true>;
template struct ausaxs::symmetry::detail::SymmetryData<true>;
