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
        std::cout << "ATOMS CHECKPOINT" << std::endl;
        std::list<Atom> im_atoms = image.generate_atoms(cutoff);
        std::cout << "ATOMS CHECKPOINT" << std::endl;
        atoms.splice(atoms.end(), im_atoms); // move im_atoms to end of atoms
        std::cout << "ATOMS CHECKPOINT" << std::endl;
    }

    // convert list to vector
    return culler->cull(atoms);
}

std::unique_ptr<Protein> em::PartialHistogramManager::generate_protein(double cutoff) const {
    vector<Atom> atoms = generate_atoms(cutoff);

    // sort vector so we can slice it into levels of charge density
    auto comparator = [] (const Atom& atom1, const Atom& atom2) {return atom1.occupancy < atom2.occupancy;};
    std::sort(atoms.begin(), atoms.end(), comparator);

    unsigned int charge_index = 0, atom_index = 0, current_index = 0;
    double charge = charge_levels[charge_index]; // initialize charge

    vector<Body> bodies(charge_levels.size());
    vector<Atom> current_atoms(atoms.size());
    
    while (atoms[atom_index].occupancy < cutoff) {atom_index++;} // search for first atom with charge larger than the cutoff
    while (charge < cutoff) {charge = charge_levels[++charge_index];} // search for first charge level larger than the cutoff 
    for (; atom_index < atoms.size(); atom_index++) {
        if (atoms[atom_index].occupancy < charge) {
            current_atoms[current_index++] = atoms[atom_index];
        } else {
            // create the body for this charge bin
            current_atoms.resize(current_index);
            bodies[charge_index] = Body(current_atoms);

            // prepare the next body
            current_index = 1;
            current_atoms.resize(atoms.size());
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

ScatteringHistogram em::PartialHistogramManager::get_histogram_slow(double cutoff) const {
    Protein protein(generate_atoms(cutoff));
    std::cout << "ATOMS IN SLOW APPROACH" << std::endl;
    for (const auto& a : protein.bodies[0].protein_atoms) {std::cout << a.as_pdb() << std::endl;}
    return protein.get_histogram();
}

void em::PartialHistogramManager::update_protein(double cutoff) {
    if (protein == nullptr || protein->bodies.empty()) {
        protein = generate_protein(cutoff); 
        protein->bind_body_signallers();

        unsigned int i = 0;
        for (const auto& b : protein->bodies) {
            b.changed_state();
        }

        previous_cutoff = cutoff;
        return;
    }
    std::unique_ptr<Protein> new_protein = generate_protein(cutoff);

    std::cout << "\nFAST APPROACH NEW PROTEIN" << std::endl;
    unsigned int i = 0;
    for (const auto& b : new_protein->bodies) {
        std::cout << "BODY " << i++ << std::endl;
        b.changed_state();
        for (const auto& a : b.protein_atoms) {
            std::cout << "\t" << a.as_pdb() << std::endl;
        }
    }

    if (cutoff < previous_cutoff) {
        std::cout << "\nCutoff smaller than previous cutoff. " << std::endl;
        for (unsigned int i = 0; i < charge_levels.size()-1; i++) {
            std::cout << "Charge level " << charge_levels[i] << " compared with previous cutoff " << previous_cutoff << std::endl;
            if (charge_levels[i] < previous_cutoff) {
                std::cout << "\tReplacing body " << i << std::endl;
                protein->bodies[i] = new_protein->bodies[i];
                protein->bodies[i].changed_state();
            } else {
                std::cout << "\tReplacing body " << i << std::endl;
                protein->bodies[i] = new_protein->bodies[i];
                protein->bodies[i].changed_state();
                break;
            }
        }        
    } else {
        std::cout << "\nCutoff larger than previous cutoff. " << std::endl;
        for (unsigned int i = 0; i < charge_levels.size(); i++) {
            std::cout << "Charge level " << charge_levels[i] << " compared with current cutoff " << cutoff << std::endl;
            if (cutoff < charge_levels[i]) {
                std::cout << "\tReplacing body " << i << std::endl;
                protein->bodies[i] = new_protein->bodies[i];
                protein->bodies[i].changed_state();
            }
        }
    }

    std::cout << "\nFAST APPROACH FINAL" << std::endl;
    i = 0;
    for (const auto& b : protein->bodies) {
        std::cout << "BODY " << i++ << (std::dynamic_pointer_cast<StateManager::UnboundSignaller>(b.signal) == nullptr ? " IS UNBOUND." : " IS BOUND.") << std::endl;
        b.changed_state();
        for (const auto& a : b.protein_atoms) {
            std::cout << "\t" << a.as_pdb() << std::endl;
        }
    }
}

std::shared_ptr<Protein> em::PartialHistogramManager::get_protein() const {return protein;}

void em::PartialHistogramManager::set_cutoff_levels(std::vector<double> levels) {
    charge_levels = levels;
    if (protein != nullptr) {protein->bodies.clear();} // we must reset the bodies to ensure they remain in sync with the new levels
}