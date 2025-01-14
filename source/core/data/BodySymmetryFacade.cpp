#include <data/BodySymmetryFacade.h>
#include <data/state/Signaller.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <io/Writer.h>

using namespace ausaxs::data::detail;

template<typename BODY, bool CONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::add(const detail::Symmetry& symmetry) requires (!CONST) {
    body->symmetries.emplace_back(symmetry);
    body->get_signaller()->internal_change();
    body->get_signaller()->external_change();
}

template<typename BODY, bool CONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::add(detail::Symmetry&& symmetry) requires (!CONST) {
    body->symmetries.emplace_back(std::move(symmetry));
    body->get_signaller()->internal_change();
    body->get_signaller()->external_change();
}

template<typename BODY, bool CONST>
std::vector<ausaxs::data::detail::Symmetry>& ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::get() requires (!CONST) {
    return body->symmetries;
}

template<typename BODY, bool CONST>
const std::vector<ausaxs::data::detail::Symmetry>& ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::get() const {
    return body->symmetries;
}

template<typename BODY, bool CONST>
ausaxs::data::detail::Symmetry& ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::get(unsigned int index) requires (!CONST) {
    return body->symmetries[index];
}

template<typename BODY, bool CONST>
const ausaxs::data::detail::Symmetry& ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::get(unsigned int index) const {
    return body->symmetries[index];
}

template<typename BODY, bool CONST>
ausaxs::data::Body ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::get_explicit_structure() const {
    std::vector<AtomFF> atoms = body->get_atoms();
    std::vector<Water> waters = body->size_water() == 0 ? std::vector<Water>{} : body->get_waters();
    atoms.reserve(body->size_symmetry_total()*body->size_atom());
    waters.reserve(body->size_symmetry_total()*body->size_water());

    for (const auto& symmetry : body->symmetries) {
        for (int i = 0; i < symmetry.repeat; ++i) {
            auto t = symmetry.get_transform(body->get_cm(), i+1);
            for (const auto& a : body->get_atoms()) {
                atoms.emplace_back(t(a.coordinates()), a.form_factor_type());
            }

            if (body->size_water() == 0) {continue;}
            for (const auto& w : body->get_waters()) {
                waters.emplace_back(t(w.coordinates()));
            }
        }
    }

    return Body(atoms, waters);
}

template<typename BODY, bool CONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::save(const io::File& path) const {
    auto body = get_explicit_structure();
    io::Writer::write(io::pdb::PDBStructure(body), path);
}

template class ausaxs::data::detail::BodySymmetryFacade<ausaxs::data::Body>;
template class ausaxs::data::detail::BodySymmetryFacade<const ausaxs::data::Body>;