/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/symmetry/BodySymmetryFacade.h>
#include <data/state/Signaller.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <io/Writer.h>

#include <span>

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
    std::vector<Water> waters = [this] () {
        auto w = body->get_waters(); return w.has_value() ? w.value().get() : std::vector<Water>{};
    }();

    if (body->size_symmetry() == 0) {
        return data::detail::SimpleBody(std::move(atoms), std::move(waters));
    }

    atoms.reserve((1+body->size_symmetry_total())*body->size_atom());
    waters.reserve((1+body->size_symmetry_total())*body->size_water());
    auto cm = body->get_cm();

    // static spans for iteration
    std::span<AtomFF> atom_span(atoms);
    std::span<Water> water_span(waters);
    for (const auto& symmetry : get()) {
        assert(atom_span.data() == atoms.data() && "atoms span has been reallocated and invalidated atom_span");
        assert(water_span.data() == waters.data() && "waters span has been reallocated and invalidated water_span");
        for (int i = 0; i < symmetry.repetitions; ++i) {
            auto t = symmetry.template get_transform<double>(cm, i+1);
            for (const auto& a : atom_span) {
                atoms.emplace_back(t(a.coordinates()), a.form_factor_type());
            }

            for (const auto& w : water_span) {
                waters.emplace_back(t(w.coordinates()));
            }
        }
    }
    assert(atoms.capacity() == atoms.size() && "atomic loop was not executed the expected number of times");
    assert(waters.capacity() == waters.size() && "water loop was not executed the expected number of times");

    return data::detail::SimpleBody(atoms, waters);
}

template<typename BODY, bool NONCONST>
void symmetry::detail::BodySymmetryFacade<BODY, NONCONST>::save(const io::File& path) const {
    auto body = explicit_structure();
    io::Writer::write(io::pdb::PDBStructure(Body(std::move(body.atoms), std::move(body.waters))), path);
}

template class symmetry::detail::BodySymmetryFacade<data::Body>;
template class symmetry::detail::BodySymmetryFacade<const data::Body>;