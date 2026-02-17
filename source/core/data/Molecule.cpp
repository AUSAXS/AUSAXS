// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/IHistogramManager.h>
#include <hist/histogram_manager/IPartialHistogramManager.h>
#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/state/BoundSignaller.h>
#include <data/state/UnboundSignaller.h>
#include <dataset/SimpleDataset.h>
#include <io/ExistingFile.h>
#include <constants/Constants.h>
#include <grid/Grid.h>
#include <grid/exv/ExvVolume.h>
#include <hydrate/ExplicitHydration.h>
#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/generation/GridBasedHydration.h>
#include <io/Writer.h>
#include <utility/Console.h>
#include <settings/All.h>

#include <numeric>
#include <cassert>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::data;

Molecule::Molecule() : bodies(), grid(nullptr), phm(nullptr), hydration_strategy(nullptr) {}

Molecule::Molecule(Molecule&& other) {*this = std::move(other);}

Molecule::~Molecule() = default;

Molecule::Molecule(std::vector<Body>&& bodies) : bodies(std::move(bodies)), grid(nullptr), phm(nullptr), hydration_strategy(nullptr) {
    initialize();
}

Molecule::Molecule(const std::vector<Body>& bodies) : bodies(bodies), grid(nullptr), phm(nullptr), hydration_strategy(nullptr) {
    initialize();
}

Molecule::Molecule(const io::File& input) : Molecule() {
    bodies = {Body(input)};
    initialize();
}

Molecule::Molecule(const std::vector<std::string>& input) : Molecule()  {
    for (const std::string& str : input) {
        bodies.emplace_back(str);
    }
    initialize();
}

#include <utility/Logging.h>
Molecule& Molecule::operator=(Molecule&& other) {
    logging::log("Molecule move-assign");
    if (this == &other) {return *this;}
    bodies = std::move(other.bodies);
    grid = std::move(other.grid);
    initialize(); // reinitialize since some of the members contains pointers to the old object
    return *this;
}

void Molecule::reset_histogram_manager() {
    set_histogram_manager(hist::factory::construct_histogram_manager(this));
}

void Molecule::initialize() {
    logging::log("initializing Molecule");
    if (!settings::flags::init_histogram_manager) {return;}
    reset_histogram_manager();
}

void Molecule::translate(const Vector3<double>& v) {
    for (auto& body : bodies) {
        body.translate(v);
    }
}

SimpleDataset Molecule::simulate_dataset(bool add_noise) const {
    SimpleDataset data = get_histogram()->debye_transform();
    data.reduce(settings::fit::N, true);
    data.simulate_errors();
    if (add_noise) {data.simulate_noise();}
    return data;
}

void Molecule::save(const io::File& path) const {
    io::Writer::write({*this}, path);
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

double Molecule::get_Rg(bool include_waters) const {
    Vector3<double> cm = get_cm(include_waters);
    double Rg = 0;

    // Rg is defined as the RMS average distance of each _electron_ from the center of mass, so multiply each atom by its effective charge
    for (auto& body : get_bodies()) {
        for (auto& a : body.get_atoms()) {
            Rg += cm.distance2(a.coordinates())*a.weight();
        }
    }

    // return the RMS
    double Z = get_total_atomic_charge();
    assert(Z != 0 && "Molecule::get_Rg: Division by zero. The molecule has no atoms.");
    return std::sqrt(Rg/get_total_atomic_charge());
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

double Molecule::get_volume_vdw() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.get_volume_vdw();});
}

double Molecule::get_volume_exv(double d) const {
    return grid::exv::get_volume_exv(this, d);
}

observer_ptr<grid::Grid> Molecule::create_grid() const {
    grid = std::make_unique<grid::Grid>(bodies); 
    return grid.get();
}

std::vector<AtomFF> Molecule::get_atoms() const {
    std::size_t N = std::accumulate(bodies.begin(), bodies.end(), std::size_t{0}, [] (std::size_t sum, const Body& body) {return sum + body.size_atom();});
    std::vector<AtomFF> atoms(N);
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

Vector3<double> Molecule::get_cm(bool include_water) const {
    Vector3<double> cm{0, 0, 0};
    double M = 0; // total mass

    // iterate through all constituent bodies
    for (const auto& body : bodies) {
        // iterate through their molecule atoms
        std::for_each(body.get_atoms().begin(), body.get_atoms().end(), [&M, &cm] (const auto& atom) {
            double m = constants::mass::get_mass(atom.form_factor_type());
            M += m;
            cm += atom.coordinates()*m;
        });
        if (!include_water) {continue;}

        // iterate through their hydration atoms
        auto w = body.get_waters();
        if (!w.has_value()) {continue;}
        std::for_each(w.value().get().begin(), w.value().get().end(), [&M, &cm] (const auto& water) {
            double m = constants::mass::get_mass(water.form_factor_type());
            M += m;
            cm += water.coords*m;
        });
    }
    assert(M != 0 && "Molecule::get_cm: Division by zero. The molecule has no atoms.");
    return cm/M;
}

std::vector<Water> Molecule::get_waters() const {
    std::vector<Water> waters(size_water());
    int n = 0; // current index
    for (const auto& body : bodies) {
        auto w = body.get_waters();
        if (!w.has_value()) {continue;}
        for (const auto& a : w.value().get()) {
            waters[n] = a;
            n++;
        }
    }
    return waters;
}

void Molecule::generate_new_hydration() {
    clear_hydration();
    if (hydration_strategy == nullptr) {
        hydration_strategy = hydrate::factory::construct_hydration_generator(this);
    }
    hydration_strategy->hydrate();
    signal_modified_hydration_layer();
}

observer_ptr<hydrate::HydrationStrategy> Molecule::get_hydration_generator() const {
    return hydration_strategy.get();
}

void Molecule::set_hydration_generator(std::unique_ptr<hydrate::HydrationStrategy> manager) {
    hydration_strategy = std::move(manager);
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
    for (auto& body : bodies) {
        body.clear_hydration();
    }
    if (grid != nullptr) {grid->clear_waters();}
    signal_modified_hydration_layer();
}

std::size_t Molecule::size_body() const {
    return bodies.size();
}

std::size_t Molecule::size_atom() const {
    return std::accumulate(bodies.begin(), bodies.end(), std::size_t{ 0 }, [] (std::size_t sum, const Body& body) {return sum + body.size_atom(); });
}

std::size_t Molecule::size_water() const {
    return std::accumulate(bodies.begin(), bodies.end(), std::size_t{ 0 }, [] (std::size_t sum, const Body& body) {return sum + body.size_water(); });
}

void Molecule::center() {
    Vector3<double> cm = get_cm();
    assert(cm.magnitude() != 0 && "Center of mass is zero. This is probably unintentional.");
    translate(-get_cm());
}

void Molecule::signal_modified_hydration_layer() const {
    if (phm == nullptr) {return;}

    // send signal to the histogram manager if relevant
    if (auto cast = dynamic_cast<IPartialHistogramManager*>(phm.get()); cast) {
        cast->signal_modified_hydration_layer();
    }
}

void Molecule::bind_body_signallers() {
    if (phm == nullptr) {return;}

    auto cast = dynamic_cast<hist::IPartialHistogramManager*>(phm.get());
    if (!cast) {
        // The caller requested the body signalling objects to be (re)bound, but the histogram manager
        // does not support this. To avoid leaving the bodies in a potentially dangerous state, we
        // register a new dummy signaller to all bodies. 
        for (unsigned int i = 0; i < bodies.size(); i++) {
            bodies[i].register_probe(std::make_shared<signaller::UnboundSignaller>());
        }
        return;
    }

    assert(cast->body_size == size_body() && "Molecule::bind_body_signallers: body size mismatch.");
    for (unsigned int i = 0; i < bodies.size(); i++) {
        bodies[i].register_probe(cast->get_probe(i));
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

void Molecule::set_histogram_manager(settings::hist::HistogramManagerChoice choice) {
    phm = hist::factory::construct_histogram_manager(this, choice);
    bind_body_signallers();
}

Body& Molecule::get_body(unsigned int index) {return bodies[index];}
const Body& Molecule::get_body(unsigned int index) const {return bodies[index];}

std::vector<Body>& Molecule::get_bodies() {return bodies;}

const std::vector<Body>& Molecule::get_bodies() const {return bodies;}

symmetry::detail::MoleculeSymmetryFacade Molecule::symmetry() const {
    return symmetry::detail::MoleculeSymmetryFacade(this);
}

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