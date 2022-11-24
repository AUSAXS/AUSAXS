#include <hist/detail/CompactCoordinates.h>

using namespace hist::detail;

CompactCoordinates::CompactCoordinates(const Body& body) : size(body.atoms().size()), data(size) {
    for (unsigned int i = 0; i < size; i++) {
        const Atom& a = body.atoms(i); 
        data[i] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
    }
}

CompactCoordinates::CompactCoordinates(const std::vector<Water>& atoms) : size(atoms.size()), data(size) {
    for (unsigned int i = 0; i < size; i++) {
        const Water& a = atoms[i]; 
        data[i] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
    }
}