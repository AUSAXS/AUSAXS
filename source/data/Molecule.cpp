#include <data/Molecule.h>
#include <data/record/Water.h>
#include <data/record/Atom.h>
#include <data/Body.h>
#include <data/state/BoundSignaller.h>
#include <data/state/UnboundSignaller.h>
#include <data/detail/AtomCollection.h>
#include <io/ExistingFile.h>
#include <hist/Histogram.h>
#include <fitter/HydrationFitter.h>
#include <settings/MoleculeSettings.h>
#include <settings/FitSettings.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <hist/distance_calculator/HistogramManagerFactory.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/detail/CompactCoordinates.h>
#include <constants/Constants.h>
#include <hydrate/Grid.h>
#include <hydrate/GridMember.h>
#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/culling/CullingStrategy.h>
#include <Symbols.h>

#include <numeric>

using namespace hist;
using namespace data;
using namespace data::record;

Molecule::Molecule() noexcept = default;

Molecule::Molecule(std::vector<Body>&& bodies) : bodies(std::move(bodies)) {
    initialize();
}

Molecule::Molecule(const std::vector<Body>& bodies) : Molecule(bodies, std::vector<Water>()) {}
Molecule::Molecule(const std::vector<Body>& bodies, const std::vector<Water>& hydration_atoms) : hydration_atoms(hydration_atoms), bodies(bodies) {
    initialize();
}

Molecule::Molecule(const std::vector<Atom>& molecule_atoms) : Molecule(molecule_atoms, std::vector<Water>()) {}
Molecule::Molecule(const std::vector<Atom>& molecule_atoms, const std::vector<Water>& hydration_atoms) : hydration_atoms(hydration_atoms) {
    bodies = {Body(molecule_atoms, this->hydration_atoms)};
    initialize();
}

Molecule::Molecule(const Molecule& molecule) : hydration_atoms(molecule.hydration_atoms), bodies(molecule.bodies), updated_charge(molecule.updated_charge), centered(molecule.centered) {
    initialize();
}

Molecule::Molecule(const io::File& input) {
    Body b1(input);
    bodies = {b1};
    this->get_waters() = std::move(bodies[0].get_waters());
    bodies[0].get_waters().clear();
    initialize();
}

Molecule::Molecule(const std::vector<std::string>& input) {
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

Molecule::~Molecule() = default;

void Molecule::initialize() {
    set_histogram_manager(hist::factory::construct_histogram_manager(this, settings::hist::weighted_bins));
    if (!centered && settings::molecule::center) {center();} // Centering *must* happen before generating the grid in 'update_effective_charge'!
    if (!updated_charge && settings::molecule::use_effective_charge) {update_effective_charge();}
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
    data.save(settings::general::output + "temp1.dat");
    data.reduce(settings::fit::N, true);
    data.save(settings::general::output + "temp2.dat");
    data.simulate_errors();
    data.save(settings::general::output + "temp3.dat");
    if (add_noise) {data.simulate_noise();}
    data.save(settings::general::output + "temp4.dat");
    return data;
}

void Molecule::save(const io::File& path) {
    // if there's only a single body, just save that instead
    if (bodies.size() == 1) {
        bodies[0].get_waters() = get_waters();
        bodies[0].save(path);
        return;
    }

    // otherwise we'll have to create a new file
    detail::AtomCollection file(get_atoms(), get_waters());
    file.write(path);
}

double Molecule::molar_mass() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.molar_mass();});
}

double Molecule::absolute_mass() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.absolute_mass();});
}

double Molecule::total_atomic_charge() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.total_atomic_charge();});
}

double Molecule::total_effective_charge() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.total_effective_charge();});
}

double Molecule::get_relative_charge() const {
    double V = get_volume_grid();
    double Z_molecule = total_atomic_charge();
    double Z_water = constants::charge::density::water*V;
    return Z_molecule - Z_water;
}

double Molecule::get_relative_charge_density() const {
    double V = get_volume_grid();
    double Z_molecule = total_atomic_charge();
    double Z_water = constants::charge::density::water*V;
    return (Z_molecule - Z_water)/V;
}

double Molecule::get_relative_mass_density() const {
    double V = get_volume_grid();
    double m_molecule = absolute_mass();
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
        for (const auto& a : body.get_atoms()) {
            double m = a.get_mass();
            M += m;
            cm += a.coords*m;
        }

        // iterate through their hydration atoms
        for (const auto& a : body.get_waters()) {
            double m = a.get_mass();
            M += m;
            cm += a.coords*m;
        }
    }

    // iterate through any generated hydration atoms
    for (const auto& a : get_waters()) {
        double m = a.get_mass();
        M += m;
        cm += a.coords*m;
    }

    return cm/M;
}

std::vector<Water>& Molecule::get_waters() {return hydration_atoms;}

const std::vector<Water>& Molecule::get_waters() const {return hydration_atoms;}

void Molecule::generate_new_hydration() {
    // delete the old hydration layer
    get_waters().clear();
    signal_modified_hydration_layer();

    // create the grid and hydrate it
    if (grid == nullptr) {create_grid();}
    else {grid->clear_waters();}
    get_waters() = grid->hydrate();
}

std::unique_ptr<hist::ICompositeDistanceHistogram> Molecule::get_histogram() const {
    return phm->calculate_all();
}

std::unique_ptr<hist::DistanceHistogram> Molecule::get_total_histogram() const {
    return phm->calculate();
}

observer_ptr<grid::Grid> Molecule::get_grid() const {
    return grid == nullptr ? create_grid() : grid.get();
}

void Molecule::set_grid(const grid::Grid& grid) {
    this->grid = std::make_unique<grid::Grid>(grid);
}

void Molecule::clear_grid() {
    grid = nullptr;
}

void Molecule::clear_hydration() {
    if (grid != nullptr) {grid->clear_waters();} // also clear the waters from the grid
    hydration_atoms.clear();
    signal_modified_hydration_layer();
}

std::size_t Molecule::body_size() const {
    return bodies.size();
}

std::size_t Molecule::atom_size() const {
    return std::accumulate(bodies.begin(), bodies.end(), std::size_t{ 0 }, [] (std::size_t sum, const Body& body) {return sum + body.get_atoms().size(); });
}

std::size_t Molecule::water_size() const {
    return hydration_atoms.size();
}

Water& Molecule::get_waters(unsigned int i) {return hydration_atoms[i];}

const Water& Molecule::get_water(unsigned int i) const {return hydration_atoms[i];}

std::vector<double> Molecule::debye_transform() const {
    auto data = hist::detail::CompactCoordinates(get_bodies());
    auto& q_axis = constants::axes::q_vals;

    auto contribution = [] (double qr, float w) -> double {
        if (qr < 1e-9) {
            return w;
        } else {
            return w*std::sin(qr)/qr;
        }
    };

    std::vector<double> I;
    I.reserve(q_axis.size());
    for (const auto& q : q_axis) {
        double sum = 0;
        for (unsigned int i = 0; i < data.size(); ++i) {
            unsigned int j = i+1;
            for (; j+7 < data.size(); j+=8) {
                auto res = data[i].evaluate(data[j], data[j+1], data[j+2], data[j+3], data[j+4], data[j+5], data[j+6], data[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {
                    sum += contribution(q*res.distances[k], 2*res.weights[k]);
                }
            }

            for (; j+3 < data.size(); j+=4) {
                auto res = data[i].evaluate(data[j], data[j+1], data[j+2], data[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {
                    sum += contribution(q*res.distances[k], 2*res.weights[k]);
                }
            }

            for (; j < data.size(); ++j) {
                auto res = data[i].evaluate(data[j]);
                    sum += contribution(q*res.distance, 2*res.weight);
            }

            sum += std::pow(data[i].value.w, 2);
        }

        sum *= std::exp(-q*q);
        I.push_back(sum);
    }
    return I;
}

void Molecule::update_effective_charge(double scaling) {
    static double previous_charge = 0;

    double displaced_vol = scaling*get_volume_grid();
    double displaced_charge = constants::charge::density::water*displaced_vol - previous_charge;
    previous_charge += displaced_charge;

    // number of atoms
    std::size_t N = atom_size();
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
    #if DEBUG
        if (phm == nullptr) [[unlikely]] {throw except::nullptr_error("Molecule::signal_modified_hydration_layer: Somehow the histogram manager has not been initialized.");}
    #endif
    phm->signal_modified_hydration_layer();
}

void Molecule::bind_body_signallers() {
    #if DEBUG
        if (phm == nullptr) [[unlikely]] {throw except::nullptr_error("Molecule::bind_body_signallers: Somehow the histogram manager has not been initialized.");}
    #endif
    for (unsigned int i = 0; i < bodies.size(); i++) {
        bodies[i].register_probe(phm->get_probe(i));
    }
}

std::shared_ptr<fitter::Fit> Molecule::fit(const io::ExistingFile& measurement) const {
    fitter::HydrationFitter fitter(measurement, get_histogram());
    return fitter.fit();
}

observer_ptr<IHistogramManager> Molecule::get_histogram_manager() const {return phm.get();}

void Molecule::set_histogram_manager(std::unique_ptr<hist::IHistogramManager> manager) {
    phm = std::move(manager);
    bind_body_signallers();
}

void Molecule::generate_unit_cell() {
    if (grid == nullptr) {create_grid();}
        // auto[min, max] = grid->bounding_box();
        Vector3<double> min, max;

    // expand box by 10%
    for (auto& v : min) {
        if (v < 0) {v *= (1 + settings::grid::scaling);} // if v is smaller than 0, multiply by 1+s
        else {      v *= (1 - settings::grid::scaling);} //                    else multiply by 1-s
    }
    for (auto& v : max) {
        if (v > 0) {v *= (1 + settings::grid::scaling);} // if v is larger than 0, multiply by 1+s
        else {      v *= (1 - settings::grid::scaling);} //                   else multiply by 1-s
    }
    auto cell_w = max - min;
    translate(-min);

    // create unit cell
    auto& file = bodies[0].get_file();
    file.header.remove("CRYST1");
    std::stringstream ss;
    ss  << "CRYST1"                                // 1 - 6
        << std::right << std::setw(8) << cell_w[0] // 7 - 15
        << std::right << std::setw(8) << cell_w[1] // 16 - 24
        << std::right << std::setw(8) << cell_w[2] // 25 - 33
        << std::right << std::setw(6) << "90"      // 34 - 40
        << std::right << std::setw(6) << "90"      // 41 - 47
        << std::right << std::setw(6) << "90"      // 48 - 54
        << " "
        << std::right << std::setw(10) << "1"      // 56 - 66
        << std::right << std::setw(4) << "P 1";    // 67 - 70
    file.add(RecordType::HEADER, ss.str());
}

void Molecule::remove_disconnected_atoms(double min_percent) {
    if (grid == nullptr) {create_grid();}
    int min = static_cast<int>(std::round(min_percent*atom_size()));
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

bool Molecule::operator==(const Molecule& other) const = default;

bool Molecule::equals_content(const Molecule& other) const {
    if (body_size() != body_size()) {return false;}
    for (unsigned int i = 0; i < body_size(); i++) {
        if (get_body(i).equals_content(other.get_body(i))) {return false;}
    }
    return true;
}