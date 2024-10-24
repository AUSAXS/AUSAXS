/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include "hist/distance_calculator/IHistogramManager.h"
#include "hist/intensity_calculator/ICompositeDistanceHistogram.h"
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Water.h>
#include <data/record/Atom.h>
#include <data/state/BoundSignaller.h>
#include <data/state/UnboundSignaller.h>
#include <data/detail/AtomCollection.h>
#include <dataset/SimpleDataset.h>
#include <io/ExistingFile.h>
#include <settings/MoleculeSettings.h>
#include <settings/FitSettings.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <hist/distance_calculator/HistogramManagerFactory.h>
#include <constants/Constants.h>
#include <grid/Grid.h>
#include <hydrate/ExplicitHydration.h>
#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/generation/GridBasedHydration.h>
#include <utility/Console.h>

#include <numeric>
#include <cassert>

using namespace hist;
using namespace data;
using namespace data::record;

Molecule::Molecule() : hydration(std::make_unique<hydrate::ExplicitHydration>()), bodies(), grid(nullptr), phm(nullptr), histogram(nullptr), hydration_strategy(nullptr) {}

Molecule::Molecule(Molecule&& other) {*this = std::move(other);}

Molecule::~Molecule() = default;

Molecule::Molecule(std::vector<Body>&& bodies) : hydration(std::make_unique<hydrate::ExplicitHydration>()), bodies(std::move(bodies)), grid(nullptr), 
                                                 phm(nullptr), histogram(nullptr), hydration_strategy(nullptr)
{
    initialize();
}

Molecule::Molecule(const std::vector<Body>& bodies) : Molecule(bodies, std::vector<Water>()) {}
Molecule::Molecule(const std::vector<Body>& bodies, const std::vector<Water>& hydration_atoms) : hydration(std::make_unique<hydrate::ExplicitHydration>(hydration_atoms)), bodies(bodies), 
                                                                                                 grid(nullptr), phm(nullptr), histogram(nullptr), hydration_strategy(nullptr) 
{
    initialize();
}

Molecule::Molecule(const std::vector<Atom>& molecule_atoms) : Molecule(molecule_atoms, std::vector<Water>()) {}
Molecule::Molecule(const std::vector<Atom>& molecule_atoms, const std::vector<Water>& hydration_atoms) : hydration(std::make_unique<hydrate::ExplicitHydration>(hydration_atoms)),
                                                                                                         grid(nullptr), phm(nullptr), histogram(nullptr), hydration_strategy(nullptr)
{
    bodies = {Body(molecule_atoms, get_waters())};
    initialize();
}

Molecule::Molecule(const io::File& input) : Molecule() {
    Body b1(input);
    bodies = {std::move(b1)};
    this->get_waters() = std::move(bodies[0].get_waters());
    bodies[0].get_waters().clear();
    initialize();
}

Molecule::Molecule(const std::vector<std::string>& input) : Molecule()  {
    std::vector<Water> waters;
    for (const std::string& str : input) {
        bodies.emplace_back(str);
        std::vector<Water>& bodyWaters = bodies.back().get_waters();
        waters.insert(waters.end(), bodyWaters.begin(), bodyWaters.end());
        bodyWaters.clear();
    }
    this->get_waters() = std::move(waters);
    initialize();
}

Molecule& Molecule::operator=(Molecule&& other) {
    if (this == &other) {return *this;}
    hydration = std::move(other.hydration);
    bodies = std::move(other.bodies);
    updated_charge = other.updated_charge;
    centered = other.centered;
    grid = std::move(other.grid);
    initialize(); // reinitialize since some of the members contains pointers to the old object
    return *this;
}

void Molecule::initialize() {
    set_histogram_manager(hist::factory::construct_histogram_manager(this, settings::hist::weighted_bins));
    if (!centered && settings::molecule::center) {center();} // Centering *must* happen before generating the grid in 'update_effective_charge'!
    if (!updated_charge && settings::molecule::use_effective_charge) {update_effective_charge();}
}

void Molecule::add_implicit_hydrogens() {
    // sanity check: if the structure already contains hydrogens, don't implicitly add more
    for (auto a : get_atoms()) {
        if (a.get_element() != constants::atom_t::H) {continue;}
        console::print_warning("Molecule::add_implicit_hydrogens: The molecule already contains hydrogen atoms. Skipping implicit addition.");
        return;
    }

    console::print_text("\tAdding implicit hydrogens to the molecule.");
    for (auto& body : bodies) {
        body.add_implicit_hydrogens();
    }
}

void Molecule::translate(const Vector3<double>& v) {
    for (auto& body : bodies) {
        body.translate(v);
    }
    for (auto& water : get_waters()) {
        water.translate(v);
    }
}

SimpleDataset Molecule::simulate_dataset(bool add_noise) const {
    SimpleDataset data = get_histogram()->debye_transform();
    data.reduce(settings::fit::N, true);
    data.simulate_errors();
    if (add_noise) {data.simulate_noise();}
    return data;
}

void Molecule::save(const io::File& path) {
    // if there's only a single body, just save that instead to preserve the original pdb header & footer
    if (bodies.size() == 1) {
        bodies[0].get_waters() = get_waters();
        bodies[0].save(path);
        return;
    }

    // otherwise we'll have to create a new file
    detail::AtomCollection file(get_atoms(), get_waters());
    file.write(path);
}

double Molecule::get_molar_mass() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.get_molar_mass();});
}

double Molecule::get_absolute_mass() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.get_absolute_mass();});
}

double Molecule::get_excluded_volume_mass() const {
    return get_volume_grid()*constants::SI::volume::A3*constants::mass::density::protein/constants::SI::mass::u;
}

double Molecule::get_total_atomic_charge() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.get_total_atomic_charge();});
}

double Molecule::get_total_effective_charge() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.get_total_effective_charge();});
}

double Molecule::get_Rg() const {
    Vector3<double> cm = get_cm();
    double Rg = 0;

    // Rg is defined as the RMS average distance of each _electron_ from the center of mass, so multiply each atom by its effective charge
    for (auto& body : get_bodies()) {
        for (auto& a : body.get_atoms()) {
            Rg += cm.distance2(a.get_coordinates())*a.get_effective_charge();
        }
    }

    // return the RMS
    return std::sqrt(Rg/get_total_effective_charge());
}

double Molecule::get_relative_charge() const {
    double V = get_volume_grid();
    double Z_molecule = get_total_atomic_charge();
    double Z_water = constants::charge::density::water*V;
    return Z_molecule - Z_water;
}

double Molecule::get_relative_charge_density() const {
    double V = get_volume_grid();
    double Z_molecule = get_total_atomic_charge();
    double Z_water = constants::charge::density::water*V;
    return (Z_molecule - Z_water)/V;
}

double Molecule::get_relative_mass_density() const {
    double V = get_volume_grid();
    double m_molecule = get_absolute_mass();
    double m_water = constants::mass::density::water*V;
    return (m_molecule - m_water)/V;
}

double Molecule::get_volume_grid() const {
    if (grid == nullptr) {create_grid();}
    return grid->get_volume();
}

observer_ptr<grid::Grid> Molecule::create_grid() const {
    grid = std::make_unique<grid::Grid>(bodies); 
    return grid.get();
}

std::vector<Atom> Molecule::get_atoms() const {
    std::size_t N = std::accumulate(bodies.begin(), bodies.end(), std::size_t{ 0 }, [] (std::size_t sum, const Body& body) { return sum + body.get_atoms().size(); });
    std::vector<Atom> atoms(N);
    int n = 0; // current index
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            atoms[n] = a;
            n++;
        }
    }
    if (n != static_cast<int>(N)) [[unlikely]] {throw except::size_error("Molecule::atoms: incorrect number of atoms. This should never happen.");}
    return atoms;
}

Vector3<double> Molecule::get_cm() const {
    Vector3<double> cm;
    double M = 0; // total mass

    // iterate through all constituent bodies
    for (const auto& body : bodies) {
        // iterate through their molecule atoms
        std::for_each(body.get_atoms().begin(), body.get_atoms().end(), [&M, &cm] (const auto& atom) {
            double m = atom.get_mass();
            M += m;
            cm += atom.coords*m;
        });

        // iterate through their hydration atoms
        std::for_each(body.get_waters().begin(), body.get_waters().end(), [&M, &cm] (const auto& water) {
            double m = water.get_mass();
            M += m;
            cm += water.coords*m;
        });
    }

    // iterate through any generated hydration atoms
    std::for_each(get_waters().begin(), get_waters().end(), [&M, &cm] (const auto& water) {
        double m = water.get_mass();
        M += m;
        cm += water.coords*m;
    });

    return cm/M;
}

const std::vector<Water>& Molecule::get_waters() const {
    if (auto cast = dynamic_cast<hydrate::ExplicitHydration*>(hydration.get()); cast != nullptr) {
        return cast->waters;
    } else {
        throw std::runtime_error("Molecule::get_waters: The chosen hydration strategy is not explicit, and thus does not generate dummy water atoms.");
    }
}

std::vector<Water>& Molecule::get_waters() {
    return const_cast<std::vector<Water>&>(static_cast<const Molecule*>(this)->get_waters());
}

observer_ptr<hydrate::HydrationStrategy> Molecule::get_hydration_generator() const {
    assert(hydration_strategy != nullptr && "Molecule::get_hydration_generator: hydration_strategy is nullptr.");
    return hydration_strategy.get();
}

void Molecule::set_hydration_generator(std::unique_ptr<hydrate::HydrationStrategy> manager) {
    hydration_strategy = std::move(manager);
}

void Molecule::generate_new_hydration() {
    if (hydration_strategy == nullptr) {
        hydration_strategy = hydrate::factory::construct_hydration_generator(this);
    }
    hydration = hydration_strategy->hydrate();
    phm->signal_modified_hydration_layer();
}

std::unique_ptr<hist::ICompositeDistanceHistogram> Molecule::get_histogram() const {
    assert(phm != nullptr && "Molecule::get_histogram: phm is nullptr.");
    return phm->calculate_all();
}

std::unique_ptr<hist::DistanceHistogram> Molecule::get_total_histogram() const {
    assert(phm != nullptr && "Molecule::get_total_histogram: phm is nullptr.");
    return phm->calculate();
}

observer_ptr<grid::Grid> Molecule::get_grid() const {
    return grid == nullptr ? create_grid() : grid.get();
}

void Molecule::set_grid(grid::Grid&& grid) {
    this->grid = std::make_unique<grid::Grid>(std::move(grid));
}

void Molecule::set_grid(std::unique_ptr<grid::Grid> grid) {
    this->grid = std::move(grid);
}

void Molecule::clear_grid() {
    grid = nullptr;
}

void Molecule::clear_hydration() {
    if (grid != nullptr) {grid->clear_waters();} // also clear the waters from the grid
    hydration->clear();
    signal_modified_hydration_layer();
}

std::size_t Molecule::size_body() const {
    return bodies.size();
}

std::size_t Molecule::size_atom() const {
    return std::accumulate(bodies.begin(), bodies.end(), std::size_t{ 0 }, [] (std::size_t sum, const Body& body) {return sum + body.get_atoms().size(); });
}

std::size_t Molecule::size_water() const {
    return get_waters().size();
}

Water& Molecule::get_waters(unsigned int i) {return get_waters()[i];}

const Water& Molecule::get_water(unsigned int i) const {return get_waters()[i];}

void Molecule::update_effective_charge(double scaling) {
    static double previous_charge = 0;

    std::size_t N = size_atom();
    if (N == 0) [[unlikely]] {return;}

    double displaced_vol = scaling*get_volume_grid();
    double displaced_charge = constants::charge::density::water*displaced_vol - previous_charge;
    previous_charge += displaced_charge;

    double charge_per_atom = -displaced_charge/N;
    if (settings::general::verbose) {
        std::cout << "Total displaced charge: " << displaced_charge << std::endl;
        std::cout << "Added " << charge_per_atom << " electrons to each atom (N: " << N << ")." << std::endl;
    }

    // subtract the charge from all molecule atoms
    for (auto& body : bodies) {
        body.update_effective_charge(charge_per_atom);
    }

    updated_charge = true;
}

void Molecule::center() {
    translate(-get_cm());
    centered = true;
}

void Molecule::signal_modified_hydration_layer() const {
    if (phm == nullptr) {return;}
    phm->signal_modified_hydration_layer();
}

void Molecule::bind_body_signallers() {
    if (phm == nullptr) {return;}
    for (unsigned int i = 0; i < bodies.size(); i++) {
        bodies[i].register_probe(phm->get_probe(i));
    }
}

observer_ptr<IHistogramManager> Molecule::get_histogram_manager() const {
    assert(phm != nullptr && "Molecule::get_histogram_manager: phm is nullptr.");
    return phm.get();
}

void Molecule::set_histogram_manager(std::unique_ptr<hist::IHistogramManager> manager) {
    phm = std::move(manager);
    bind_body_signallers();
}

void Molecule::remove_disconnected_atoms(double min_percent) {
    if (grid == nullptr) {create_grid();}
    int min = static_cast<int>(std::round(min_percent*size_atom()));
    auto to_remove = grid->remove_disconnected_atoms(min);

    // sanity check
    if (to_remove.size() != get_atoms().size()) {
        throw except::unexpected("Molecule::remove_disconnected_atoms: "
        "The number of atoms to remove (" + std::to_string(to_remove.size()) + ") does not match the number of molecule atoms (" + std::to_string(get_atoms().size()) + ").");
    }

    // remove the atoms from the molecule bodies
    unsigned int index = 0;
    for (auto& body : bodies) {
        unsigned int removed = 0;
        std::vector<Atom> new_atoms(body.get_atoms().size());
        std::vector<Atom>& atoms = body.get_atoms();
        for (unsigned int i = 0; i < atoms.size(); i++) {
            if (to_remove[index + i]) {
                removed++;
            } else {
                new_atoms[i-removed] = std::move(atoms[i]);
            }
        }
        new_atoms.resize(atoms.size() - removed);
        body.get_atoms() = std::move(new_atoms);
    }
}

Body& Molecule::get_body(unsigned int index) {return bodies[index];}
const Body& Molecule::get_body(unsigned int index) const {return bodies[index];}

std::vector<Body>& Molecule::get_bodies() {return bodies;}

const std::vector<Body>& Molecule::get_bodies() const {return bodies;}

bool Molecule::equals_content(const Molecule& other) const {
    if (size_body() != other.size_body()) {
        return false;
    }

    for (unsigned int i = 0; i < size_body(); i++) {
        if (get_body(i).equals_content(other.get_body(i))) {
            return false;
        }
    }
    return true;
}