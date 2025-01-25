#include <data/symmetry/BodySymmetryFacade.h>
#include <data/state/Signaller.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <io/Writer.h>

using namespace ausaxs::data::detail;

template<typename BODY, bool NONCONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::add(symmetry::Symmetry&& symmetry) requires (NONCONST) {
    assert(body->get_signaller() && "BodySymmetryFacade::add: Body signaller object not initialized.");
    body->get_signaller()->symmetry_changed();
    body->symmetries->get().emplace_back(std::move(symmetry));
}

template<typename BODY, bool NONCONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::add(symmetry::type symmetry) requires (NONCONST) {
    assert(body->get_signaller() && "BodySymmetryFacade::add: Body signaller object not initialized.");
    body->get_signaller()->symmetry_changed();
    body->symmetries->add(symmetry);
}

template<typename BODY, bool NONCONST>
std::vector<ausaxs::symmetry::Symmetry>& ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::get() requires (NONCONST) {
    return body->symmetries->get();
}

template<typename BODY, bool NONCONST>
const std::vector<ausaxs::symmetry::Symmetry>& ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::get() const {
    return body->symmetries->get();
}

template<typename BODY, bool NONCONST>
ausaxs::symmetry::Symmetry& ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::get(unsigned int index) requires (NONCONST) {
    assert(index < body->symmetries->get().size());
    return body->symmetries->get()[index];
}

template<typename BODY, bool NONCONST>
const ausaxs::symmetry::Symmetry& ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::get(unsigned int index) const {
    assert(index < body->symmetries->get().size());
    return body->symmetries->get()[index];
}

template<typename BODY, bool NONCONST>
ausaxs::observer_ptr<ausaxs::symmetry::SymmetryStorage> ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::get_obj() requires (NONCONST) {
    return body->symmetries.get();
}

template<typename BODY, bool NONCONST>
ausaxs::observer_ptr<const ausaxs::symmetry::SymmetryStorage> ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::get_obj() const {
    return body->symmetries.get();
}

template<typename BODY, bool NONCONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::set_obj(std::unique_ptr<ausaxs::symmetry::SymmetryStorage> obj) requires (NONCONST) {
    body->symmetries = std::move(obj);
}

template<typename BODY, bool NONCONST>
ausaxs::data::Body ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::get_explicit_structure() const {
    std::vector<AtomFF> atoms = body->get_atoms();
    std::vector<Water> waters = body->size_water() == 0 ? std::vector<Water>{} : body->get_waters();
    atoms.reserve(body->size_symmetry_total()*body->size_atom());
    waters.reserve(body->size_symmetry_total()*body->size_water());

    for (const auto& symmetry : get()) {
        for (int i = 0; i < symmetry.repeat; ++i) {
            auto t = symmetry.template get_transform<double>(i+1);
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

template<typename BODY, bool NONCONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, NONCONST>::save(const io::File& path) const {
    auto body = get_explicit_structure();
    io::Writer::write(io::pdb::PDBStructure(body), path);
}

template class ausaxs::data::detail::BodySymmetryFacade<ausaxs::data::Body>;
template class ausaxs::data::detail::BodySymmetryFacade<const ausaxs::data::Body>;