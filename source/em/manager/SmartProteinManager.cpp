#include <em/manager/SmartProteinManager.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/Body.h>
#include <utility/Console.h>
#include <em/detail/ImageStackBase.h>
#include <settings/EMSettings.h>
#include <em/Image.h>
#include <data/record/Water.h>
#include <Symbols.h>

#include <vector>
#include <functional>

using namespace em::managers;
using namespace data;
using namespace data::record;

std::unique_ptr<hist::ICompositeDistanceHistogram> SmartProteinManager::get_histogram(double cutoff) {
    update_protein(cutoff);
    return protein->get_histogram();
}

std::vector<Atom> SmartProteinManager::generate_atoms(double cutoff) const {
    // we use a list since we will have to append quite a few other lists to it
    std::list<Atom> atoms;
    const auto& imagestack = images->images();
    unsigned int step = settings::em::sample_frequency;
    for (unsigned int i = 0; i < imagestack.size(); i += step) {
        std::list<Atom> im_atoms = imagestack[i].generate_atoms(cutoff);
        atoms.splice(atoms.end(), im_atoms); // move im_atoms to end of atoms
    }

    // convert list to vector
    return std::vector<Atom>(std::make_move_iterator(std::begin(atoms)), std::make_move_iterator(std::end(atoms)));
}

std::unique_ptr<data::Molecule> SmartProteinManager::generate_protein(double cutoff) const {
    std::vector<Atom> atoms = generate_atoms(cutoff);
    std::vector<Body> bodies(charge_levels.size());
    std::vector<Atom> current_atoms(atoms.size());

    if (atoms.empty()) {
        console::print_warning("Warning in SmartProteinManager::generate_protein: No voxels found for cutoff \"" + std::to_string(cutoff) + "\".");
        return std::make_unique<data::Molecule>(bodies);
    }

    if (charge_levels.empty()) {
        throw except::out_of_bounds("SmartProteinManager::generate_protein: charge_levels is empty.");
    }

    std::function<bool(double, double)> compare_positive = [] (double v1, double v2) {return v1 < v2;};
    std::function<bool(double, double)> compare_negative = [] (double v1, double v2) {return v1 > v2;};
    std::function<bool(double, double)> compare_func = 0 <= cutoff ? compare_positive : compare_negative;

    // sort vector so we can slice it into levels of charge density
    auto comparator = [&compare_func] (const Atom& atom1, const Atom& atom2) {return compare_func(atom1.occupancy, atom2.occupancy);};
    std::sort(atoms.begin(), atoms.end(), comparator);

    unsigned int charge_index = 0, atom_index = 0, current_index = 0;
    double charge = charge_levels[charge_index]; // initialize charge

    while (compare_func(atoms[atom_index].occupancy, cutoff)) {atom_index++;} // search for first atom with charge larger than the cutoff
    while (compare_func(charge, cutoff)) {charge = charge_levels[++charge_index];} // search for first charge level larger than the cutoff 
    while (atom_index < atoms.size()) {
        if (compare_func(atoms[atom_index].occupancy, charge)) {
            current_atoms[current_index++] = atoms[atom_index++];
        } else {
            // create the body for this charge bin
            current_atoms.resize(current_index);
            bodies[charge_index] = Body(current_atoms);

            // prepare the next body
            current_index = 0;
            current_atoms.resize(atoms.size());

            // increment the charge level
            if (charge_index+1 == charge_levels.size()) [[unlikely]] {
                throw except::unexpected("SmartProteinManager::generate_protein: Reached end of charge levels list.");
            }
            charge = charge_levels[++charge_index];
        }
    }
    
    // create the final body of the loop
    current_atoms.resize(current_index);
    bodies[charge_index] = Body(current_atoms);

    return std::make_unique<data::Molecule>(bodies);
}

void SmartProteinManager::update_protein(double cutoff) {
    if (protein == nullptr || protein->atom_size() == 0) {
        protein = generate_protein(cutoff); 
        protein->bind_body_signallers();

        previous_cutoff = cutoff;
        return;
    }

    if (cutoff == previous_cutoff) {
        return;
    }

    // sanity check
    if (charge_levels.empty()) {
        throw except::unexpected("SmartProteinManager::update_protein: charge_levels is empty.");
    }

    std::unique_ptr<data::Molecule> new_protein = generate_protein(cutoff);
    // std::cout << "Found " << new_protein->atom_size() << " voxels with a cutoff larger than " << cutoff << std::endl;

    std::function<bool(double, double)> compare_positive = [] (double v1, double v2) {return v1 < v2;};
    std::function<bool(double, double)> compare_negative = [] (double v1, double v2) {return v1 > v2;};
    std::function<bool(double, double)> compare_func = 0 <= cutoff ? compare_positive : compare_negative;

    if (compare_func(cutoff, previous_cutoff)) {
        // since cutoff is smaller than previously, we have to change all bins in the range [cutoff, previous_cutoff]

        // skip all bins before the relevant range
        unsigned int charge_index = 0;
        double current_cutoff = charge_levels[0];
        while (compare_func(current_cutoff, cutoff) && charge_index < charge_levels.size()) {
            current_cutoff = charge_levels[++charge_index];}

        // iterate through the remaining bins, and use a break statement to stop when we leave the relevant range
        for (; charge_index < charge_levels.size()-1; charge_index++) {
            // std::cout << "Charge level " << charge_levels[i] << " compared with previous cutoff " << previous_cutoff << std::endl;
            // check if the current bin is inside the range
            if (compare_func(charge_levels[charge_index], previous_cutoff)) {
                // if so, we replace it with the new contents
                protein->get_body(charge_index) = std::move(new_protein->get_body(charge_index));
            } else {
                // if we have the same number of atoms as earlier, nothing has changed
                if (new_protein->get_body(charge_index).atom_size() == protein->get_body(charge_index).atom_size()) {
                    break;
                }

                // otherwise we replace it and stop iterating
                protein->get_body(charge_index) = std::move(new_protein->get_body(charge_index));
                break;
            }
        }
    } else {
        // since cutoff is larger than previously, we have to change all bins in the range [previous_cutoff, cutoff]

        // skip all bins before the relevant range
        unsigned int charge_index = 0;
        double current_cutoff = charge_levels[0];
        while (compare_func(current_cutoff, previous_cutoff) && charge_index < charge_levels.size()) {
            current_cutoff = charge_levels[++charge_index];}

        // iterate through the remaining bins, and use a break statement to stop when we leave the relevant range
        for (; charge_index < charge_levels.size()-1; charge_index++) {
            // std::cout << "Charge level " << charge_levels[i] << " compared with previous cutoff " << previous_cutoff << std::endl;
            // check if the current bin is inside the range
            if (compare_func(charge_levels[charge_index], cutoff)) {
                // if so, we replace it with the new contents
                protein->get_body(charge_index) = std::move(new_protein->get_body(charge_index));
            } else {
                // otherwise we replace it and stop iterating
                protein->get_body(charge_index) = std::move(new_protein->get_body(charge_index));
                break;
            }
        }

    }

    previous_cutoff = cutoff;
}

data::Molecule* SmartProteinManager::get_protein() const {
    #ifdef DEBUG
        if (protein == nullptr) {throw except::nullptr_error("SmartProteinManager::get_protein: Protein has not been initialized yet.");}
    #endif
    return protein.get();
}

data::Molecule* SmartProteinManager::get_protein(double cutoff) {
    update_protein(cutoff);
    return protein.get();
}

void SmartProteinManager::set_charge_levels(const std::vector<double>& levels) noexcept {
    ProteinManager::set_charge_levels(levels);
    protein = nullptr; // the protein must be generated anew to ensure the bodies remains in sync with the new levels
}