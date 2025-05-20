/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/manager/SmartProteinManager.h>
#include <em/detail/ImageStackBase.h>
#include <em/detail/EMGrid.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/Console.h>
#include <utility/Limit3D.h>
#include <utility/Logging.h>
#include <settings/HistogramSettings.h>
#include <settings/ExvSettings.h>
#include <settings/EMSettings.h>
#include <settings/GridSettings.h>
#include <settings/Flags.h>

#include <vector>
#include <cassert>

using namespace ausaxs;
using namespace ausaxs::em::managers;
using namespace ausaxs::data;

SmartProteinManager::~SmartProteinManager() = default;

std::unique_ptr<hist::ICompositeDistanceHistogram> SmartProteinManager::get_histogram(double cutoff) {
    return get_protein(cutoff)->get_histogram();
}

observer_ptr<data::Molecule> SmartProteinManager::get_protein(double cutoff) {
    update_protein(cutoff);

    auto ax = images->get_header()->get_axes();

    // ensure the grid is large enough to contain the entire map
    //? is this a waste of time since it is being overwritten for every iteration anyway?
    Limit3D limits(ax.x.min, ax.x.max, ax.y.min, ax.y.max, ax.z.min, ax.z.max);
    double expand_x = 0.5*limits.x.span()*settings::grid::scaling;
    double expand_y = 0.5*limits.y.span()*settings::grid::scaling;
    double expand_z = 0.5*limits.z.span()*settings::grid::scaling;
    limits.x.min -= expand_x, limits.x.max += expand_x;
    limits.y.min -= expand_y, limits.y.max += expand_y;
    limits.z.min -= expand_z, limits.z.max += expand_z;

    // hydration must be cleared to ensure the grid state will be exactly the same as the one used to generate the hydration
    // when the grid is first created, *only* the body atoms are added to the grid, with the hydration being added later
    // if we do not clear the hydration now, since it can be bound to individual bodies, the hydration will be added in-between the body atoms, 
    // potentially leading to a different grid state than the one used to generate the hydration
    protein->clear_hydration();

    protein->set_grid(std::make_unique<grid::EMGrid>(limits));
    for (auto& body : protein->get_bodies()) {
        protein->get_grid()->add(body);
    }
    if (settings::em::hydrate) {protein->generate_new_hydration();}
    return protein.get();
}

observer_ptr<const data::Molecule> SmartProteinManager::get_protein() const {
    assert(protein != nullptr && "SmartProteinManager::get_protein: Protein has not been initialized yet.");
    return protein.get();
}

void SmartProteinManager::set_charge_levels(const std::vector<double>& levels) noexcept {
    ProteinManager::set_charge_levels(levels);
    protein = nullptr; // the protein must be generated anew to ensure the bodies remains in sync with the new levels
}

std::vector<EMAtom> SmartProteinManager::generate_atoms(double cutoff) const {
    // we use a list since we will have to append quite a few other lists to it
    std::list<EMAtom> atoms;
    const auto& imagestack = images->images();
    unsigned int step = settings::em::sample_frequency;
    for (unsigned int i = 0; i < imagestack.size(); i += step) {
        std::list<EMAtom> im_atoms = imagestack[i].generate_atoms(cutoff);
        atoms.splice(atoms.end(), im_atoms); // move im_atoms to end of atoms
    }

    // convert list to vector
    return std::vector<EMAtom>(std::make_move_iterator(std::begin(atoms)), std::make_move_iterator(std::end(atoms)));
}

std::unique_ptr<data::Molecule> SmartProteinManager::generate_new_protein(double cutoff) const {
    std::vector<EMAtom> atoms = generate_atoms(cutoff);
    std::vector<Body> bodies(charge_levels.size());
    std::vector<AtomFF> current_atoms(atoms.size());

    if (atoms.empty()) {
        console::print_warning("Warning in SmartProteinManager::generate_protein: No voxels found for cutoff \"" + std::to_string(cutoff) + "\".");
        return std::make_unique<data::Molecule>(bodies);
    }

    if (charge_levels.empty()) {
        throw except::out_of_bounds("SmartProteinManager::generate_protein: charge_levels is empty.");
    }

    // sort vector so we can slice it into levels of charge density
    std::sort(atoms.begin(), atoms.end(), [] (const EMAtom& atom1, const EMAtom& atom2) {return atom1.charge_density() < atom2.charge_density();});

    unsigned int charge_index = 0, atom_index = 0, current_index = 0;
    double charge = charge_levels[charge_index]; // initialize charge

    while (atoms[atom_index].charge_density() < cutoff) {++atom_index;} // search for first atom with charge larger than the cutoff
    while (charge < cutoff) {charge = charge_levels[++charge_index];}   // search for first charge level larger than the cutoff 
    while (atom_index < atoms.size()) {
        if (atoms[atom_index].charge_density() < charge) {
            current_atoms[current_index++] = atoms[atom_index++].get_atom_ff();
        } else {
            // create the body for this charge bin
            current_atoms.resize(current_index);
            bodies[charge_index] = Body(current_atoms);

            // prepare the next body
            current_index = 0;
            current_atoms.resize(atoms.size() - atom_index);

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

    return std::make_unique<data::Molecule>(std::move(bodies));
}

void SmartProteinManager::toggle_histogram_manager_init(bool state) {
    settings::flags::init_histogram_manager = state;
}

void SmartProteinManager::update_protein(double cutoff) {
    if (protein == nullptr || protein->size_atom() == 0) {
        // the protein is not initialized, so simply assign it a new one
        toggle_histogram_manager_init(false);
        logging::log("SmartProteinManager::update_protein: protein is nullptr or empty. Generating new protein.");
        protein = generate_new_protein(cutoff); 
        protein->set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManagerMT);
        previous_cutoff = cutoff;
        return;
    }

    if (cutoff == previous_cutoff) {
        logging::log("SmartProteinManager::update_protein: cutoff is the same as previous. Nothing to do.");
        return;
    }

    // sanity check
    if (charge_levels.empty()) {
        throw except::unexpected("SmartProteinManager::update_protein: charge_levels is empty.");
    }
    logging::log("SmartProteinManager::update_protein: cutoff = " + std::to_string(cutoff) + ", previous_cutoff = " + std::to_string(previous_cutoff));

    // the idea is to generate a new protein with the new cutoff, and then replace the relevant bodies in the old protein with the new ones
    // this will trigger the "update" signals for the replaced bodies, which will inform the histogram manager to only recalculate their contributions
    std::unique_ptr<data::Molecule> new_protein = generate_new_protein(cutoff);
    if (cutoff < previous_cutoff) {
        // since cutoff is smaller than previously, we have to change all bins in the range [cutoff, previous_cutoff]

        // skip all bins before the relevant range
        unsigned int charge_index = 0;
        double current_cutoff = charge_levels[0];
        while (current_cutoff < cutoff && charge_index < charge_levels.size()) {
            current_cutoff = charge_levels[++charge_index];
        }

        // iterate through the remaining bins, and use a break statement to stop when we leave the relevant range
        for (; charge_index < charge_levels.size(); ++charge_index) {
            // check if the current bin is inside the range
            if (charge_levels[charge_index] < previous_cutoff) {
                // if so, we replace it with the new contents
                protein->get_body(charge_index) = std::move(new_protein->get_body(charge_index));
            } else {
                // if we have the same number of atoms as earlier, nothing has changed
                if (new_protein->get_body(charge_index).size_atom() == protein->get_body(charge_index).size_atom()) {
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
        while (current_cutoff < previous_cutoff && charge_index < charge_levels.size()) {
            current_cutoff = charge_levels[++charge_index];
        }

        // iterate through the remaining bins, and use a break statement to stop when we leave the relevant range
        for (; charge_index < charge_levels.size(); ++charge_index) {
            // check if the current bin is inside the range
            if (charge_levels[charge_index] < cutoff) {
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