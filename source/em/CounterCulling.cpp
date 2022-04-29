#include <em/CounterCulling.h>

using std::vector;

vector<Atom> em::CounterCulling::cull(std::list<Atom>& atoms) const {
    double percent = 1 - double(target_count)/atoms.size();

    // prepare output vector
    vector<Atom> output;
    output.reserve(atoms.size());

    // iterate through the input list
    double counter = 0;
    for (const Atom& atom : atoms) {
        counter += percent;
        if (counter > 1) {counter--; continue;}
        output.push_back(std::move(atom));
    }
    return output;
}