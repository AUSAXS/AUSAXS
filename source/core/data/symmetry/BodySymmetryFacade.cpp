#include <data/symmetry/BodySymmetryFacade.h>
#include <data/state/Signaller.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <io/Writer.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::data::detail;

template<typename BODY, bool NONCONST>
void symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::add(symmetry::Symmetry&& symmetry) requires (NONCONST) {
    assert(body->get_signaller() && "BodySymmetryFacade::add: Body signaller object not initialized.");
    body->symmetries->get().emplace_back(std::move(symmetry));
    body->get_signaller()->set_symmetry_size(body->size_symmetry());
}

template<typename BODY, bool NONCONST>
void symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::add(symmetry::type symmetry) requires (NONCONST) {
    assert(body->get_signaller() && "BodySymmetryFacade::add: Body signaller object not initialized.");
    body->symmetries->add(symmetry);
    body->get_signaller()->set_symmetry_size(body->size_symmetry());
}

template<typename BODY, bool NONCONST>
std::vector<symmetry::Symmetry>& symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::get() requires (NONCONST) {
    for (std::size_t i = 0; i < body->size_symmetry(); ++i) {body->get_signaller()->modified_symmetry(i);}
    return body->symmetries->get();
}

template<typename BODY, bool NONCONST>
const std::vector<symmetry::Symmetry>& symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::get() const {
    return body->symmetries->get();
}

template<typename BODY, bool NONCONST>
symmetry::Symmetry& symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::get(unsigned int index) requires (NONCONST) {
    assert(index < body->symmetries->get().size());
    body->get_signaller()->modified_symmetry(index);
    return body->symmetries->get()[index];
}

template<typename BODY, bool NONCONST>
const symmetry::Symmetry& symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::get(unsigned int index) const {
    assert(index < body->symmetries->get().size());
    return body->symmetries->get()[index];
}

template<typename BODY, bool NONCONST>
observer_ptr<symmetry::SymmetryStorage> symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::get_obj() requires (NONCONST) {
    for (std::size_t i = 0; i < body->size_symmetry(); ++i) {body->get_signaller()->modified_symmetry(i);}
    return body->symmetries.get();
}

template<typename BODY, bool NONCONST>
observer_ptr<const symmetry::SymmetryStorage> symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::get_obj() const {
    return body->symmetries.get();
}

template<typename BODY, bool NONCONST>
void symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::set_obj(std::unique_ptr<symmetry::SymmetryStorage> obj) requires (NONCONST) {
    body->symmetries = std::move(obj);
}

template<typename BODY, bool NONCONST>
std::size_t symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::size_atom_total() const {
    return body->size_atom()*(body->size_symmetry_total()+1);
}

template<typename BODY, bool NONCONST>
std::size_t symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::size_water_total() const {
    return body->size_water()*(body->size_symmetry_total()+1);
}

template<typename BODY, bool NONCONST>
data::detail::SimpleBody symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::explicit_structure() const {
    std::vector<AtomFF> atoms = body->get_atoms();
    std::vector<Water> waters = body->size_water() == 0 ? std::vector<Water>{} : body->get_waters();

    if (body->size_symmetry() == 0) {
        return data::detail::SimpleBody(std::move(atoms), std::move(waters));
    }

    atoms.reserve(body->size_symmetry_total()*body->size_atom());
    waters.reserve(body->size_symmetry_total()*body->size_water());
    auto cm = body->get_cm();

    for (const auto& symmetry : get()) {
        for (int i = 0; i < symmetry.repetitions; ++i) {
            auto t = symmetry.template get_transform<double>(cm, i+1);
            for (const auto& a : body->get_atoms()) {
                atoms.emplace_back(t(a.coordinates()), a.form_factor_type());
            }

            if (body->size_water() == 0) {continue;}
            for (const auto& w : body->get_waters()) {
                waters.emplace_back(t(w.coordinates()));
            }
        }
    }

    return data::detail::SimpleBody(atoms, waters);
}

template<typename BODY, bool NONCONST>
void symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::save(const io::File& path) const {
    auto body = explicit_structure();
    io::Writer::write(io::pdb::PDBStructure(Body(std::move(body.atoms), std::move(body.waters))), path);
}

template class symmetry::detail::BodySymmetryFacade<data::Body>;
template class symmetry::detail::BodySymmetryFacade<const data::Body>;