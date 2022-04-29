#include <em/NoCulling.h>

std::vector<Atom> em::NoCulling::cull(std::list<Atom>& atoms) const {
    std::vector<Atom> output;
    output.reserve(atoms.size());
    output.assign(std::move_iterator(atoms.begin()), std::move_iterator(atoms.end()));
    return output;
}