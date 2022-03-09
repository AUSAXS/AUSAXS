#include <em/NoCulling.h>

vector<Atom> em::NoCulling::cull(list<Atom>& atoms) const {
    vector<Atom> output;
    output.reserve(atoms.size());
    output.assign(std::move_iterator(atoms.begin()), std::move_iterator(atoms.end()));
    return output;
}