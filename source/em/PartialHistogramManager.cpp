#include <vector>

#include <em/PartialHistogramManager.h>
#include <em/ImageStack.h>
#include <em/CounterCulling.h>
#include <ScatteringHistogram.h>
#include <data/Protein.h>

em::PartialHistogramManager::PartialHistogramManager(const ImageStack& images) :images(images), culler(std::make_unique<CounterCulling>()) {}

ScatteringHistogram em::PartialHistogramManager::get_histogram(double cutoff) {
    update_protein(cutoff);
    return protein->get_histogram();
}

std::vector<Atom> em::PartialHistogramManager::generate_atoms(double cutoff) const {
    // we use a list since we will have to append quite a few other lists to it
    std::list<Atom> atoms;
    for (const Image& image: images.images()) {
        std::list<Atom> im_atoms = image.generate_atoms(cutoff);
        atoms.splice(atoms.end(), im_atoms); // move im_atoms to end of atoms
    }

    // convert list to vector
    return culler->cull(atoms);
}

std::unique_ptr<Protein> em::PartialHistogramManager::generate_protein(double cutoff) const {
    vector<Atom> atoms = generate_atoms(cutoff);

    // sort vector so we can slice it into levels of charge density
    auto comparator = [] (const Atom& atom1, const Atom& atom2) {return atom1.occupancy < atom2.occupancy;};
    std::sort(atoms.begin(), atoms.end(), comparator);

    unsigned int charge_index = 0, current_index = 0;
    double charge = charge_levels[charge_index]; // initialize charge

    vector<Body> bodies(charge_levels.size());
    vector<Atom> current_atoms(atoms.size()/2);
    for (unsigned int atom_index = 0; atom_index < atoms.size(); atom_index++) {
        if (atoms[atom_index].occupancy < charge) {
            current_atoms[current_index++] = atoms[atom_index];
        } else {
            // create the body for this charge bin
            current_atoms.resize(current_index);
            bodies[charge_index] = Body(current_atoms);

            // prepare the next body
            current_index = 1;
            current_atoms.resize(atoms.size()/2);
            current_atoms[0] = atoms[atom_index]; // add the atom of the current iteration

            // increment the charge level
            charge = charge_levels[++charge_index];
        }
    }
    // create the final body of the loop
    current_atoms.resize(current_index);
    bodies[charge_index] = Body(current_atoms);

    return std::make_unique<Protein>(bodies);
}

void em::PartialHistogramManager::update_protein(double cutoff) {
    if (protein == nullptr || protein->bodies.empty()) {std::cout << "Initializing primary protein." << std::endl; protein = generate_protein(cutoff); return;}
    std::unique_ptr<Protein> new_protein = generate_protein(cutoff);

    for (unsigned int i = 0; i < charge_levels.size(); i++) {
        if (charge_levels[i] < cutoff) {
            protein->bodies[i].protein_atoms.clear();
            protein->bodies[i].changed_state();
        } else if (i+1 < charge_levels.size() && cutoff < charge_levels[i+1]) {
            protein->bodies[i] = new_protein->bodies[i];
            protein->bodies[i].changed_state();
            break;
        }
    }
}

std::shared_ptr<Protein> em::PartialHistogramManager::get_protein() const {return protein;}

void em::PartialHistogramManager::set_cutoff_levels(std::vector<double> levels) {
    charge_levels = levels;
    if (protein != nullptr) {protein->bodies.clear();} // we must reset the bodies to ensure they remain in sync with the new levels
}