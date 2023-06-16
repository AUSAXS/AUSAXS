#include <data/Protein.h>
#include <data/Water.h>
#include <data/Atom.h>
#include <data/Body.h>
#include <data/state/BoundSignaller.h>
#include <data/state/UnboundSignaller.h>
#include <io/ProteinFile.h>
#include <io/ExistingFile.h>
#include <hist/Histogram.h>
#include <fitter/HydrationFitter.h>
#include <settings/ProteinSettings.h>
#include <settings/FitSettings.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <hist/detail/HistogramManagerFactory.h>
#include <hist/HistogramManager.h>
#include <utility/Constants.h>
#include <hydrate/Grid.h>
#include <hydrate/GridMember.h>
#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/culling/CullingStrategy.h>

using namespace hist;

Protein::Protein(std::vector<Body>&& bodies) : bodies(std::move(bodies)) {
    initialize();
}

Protein::Protein(const std::vector<Body>& bodies) : Protein(bodies, std::vector<Water>()) {}
Protein::Protein(const std::vector<Body>& bodies, const std::vector<Water>& hydration_atoms) : hydration_atoms(hydration_atoms), bodies(bodies) {
    initialize();
}

Protein::Protein(const std::vector<Atom>& protein_atoms) : Protein(protein_atoms, std::vector<Water>()) {}
Protein::Protein(const std::vector<Atom>& protein_atoms, const std::vector<Water>& hydration_atoms) : hydration_atoms(hydration_atoms) {
    bodies = {Body(protein_atoms, this->hydration_atoms)}; // 'this' keyword is necessary, otherwise the objects are bound to the argument instead of the member
    initialize();
}

Protein::Protein(const std::vector<std::vector<Atom>>& protein_atoms) : Protein(protein_atoms, std::vector<Water>()) {}
Protein::Protein(const std::vector<std::vector<Atom>>& protein_atoms, const std::vector<Water>& hydration_atoms) : hydration_atoms(hydration_atoms) {
    for (unsigned int i = 0; i < protein_atoms.size(); i++) {
        bodies.push_back(Body(protein_atoms[i], std::vector<Water>(0)));
    }
    initialize();
}

Protein::Protein(const Protein& protein) : hydration_atoms(protein.hydration_atoms), bodies(protein.bodies), updated_charge(protein.updated_charge), centered(protein.centered) {
    initialize();
}

Protein::Protein(const io::ExistingFile& input) {
    Body b1(input);
    bodies = {b1};
    this->get_waters() = std::move(bodies[0].get_waters());
    bodies[0].get_waters().clear();
    initialize();
}

Protein::Protein(const std::vector<std::string>& input) {
    for (size_t i = 0; i < input.size(); i++) {
        bodies.push_back(Body(input[i]));
    }
    initialize();
}

Protein::~Protein() = default;

void Protein::initialize() {
    set_histogram_manager(hist::factory::construct_histogram_manager(this));
    create_grid();
    if (settings::protein::use_effective_charge) {
        update_effective_charge();
    }
}

void Protein::translate(const Vector3<double>& v) {
    for (auto& body : bodies) {
        body.translate(v);
    }
    for (auto& hetatom : get_waters()) {
        hetatom.translate(v);
    }
}

SimpleDataset Protein::simulate_dataset(bool add_noise) const {
    SimpleDataset data = get_histogram().calc_debye_scattering_intensity();
    data.reduce(settings::fit::N, true);
    data.simulate_errors();
    if (add_noise) {
        data.simulate_noise();
    }
    return data;
}

void Protein::save(const io::File& path) {
    // if there's only a single body, just save that instead
    if (bodies.size() == 1) {
        bodies[0].get_waters() = get_waters();
        bodies[0].save(path);
        return;
    }

    // otherwise we'll have to create a new file
    ProteinFile file(get_atoms(), get_waters());
    file.write(path);
}

double Protein::get_volume_acids() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.get_volume_acids();});
}

double Protein::molar_mass() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.molar_mass();});
}

double Protein::absolute_mass() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.absolute_mass();});
}

double Protein::total_atomic_charge() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.total_atomic_charge();});
}

double Protein::total_effective_charge() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.total_effective_charge();});
}

double Protein::get_relative_charge() const {
    double V = get_volume_grid();
    double Z_protein = total_atomic_charge();
    double Z_water = constants::charge::density::water*V;
    return Z_protein - Z_water;
}

double Protein::get_relative_charge_density() const {
    double V = get_volume_grid();
    double Z_protein = total_atomic_charge();
    double Z_water = constants::charge::density::water*V;
    return (Z_protein - Z_water)/V;
}

double Protein::get_relative_mass_density() const {
    double V = get_volume_grid();
    double m_protein = absolute_mass();
    double m_water = constants::mass::density::water*V;
    return (m_protein - m_water)/V;
}

double Protein::get_volume_grid() const {
    return grid->get_volume();
}

std::shared_ptr<grid::Grid> Protein::create_grid() {
    grid = std::make_shared<grid::Grid>(bodies); 
    return grid;
}

std::vector<Atom> Protein::get_atoms() const {
    int N = std::accumulate(bodies.begin(), bodies.end(), 0, [] (double sum, const Body& body) {return sum + body.get_atoms().size();});
    std::vector<Atom> atoms(N);
    int n = 0; // current index
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            atoms[n] = a;
            n++;
        }
    }
    if (n != N) [[unlikely]] {throw except::size_error("Protein::atoms: incorrect number of atoms. This should never happen.");}
    return atoms;
}

Vector3<double> Protein::get_cm() const {
    Vector3<double> cm;
    double M = 0; // total mass

    // iterate through all constituent bodies
    for (const auto& body : bodies) {
        // iterate through their protein atoms
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

std::vector<Water>& Protein::get_waters() {return hydration_atoms;}

const std::vector<Water>& Protein::get_waters() const {return hydration_atoms;}

void Protein::generate_new_hydration() {
    // delete the old hydration layer
    get_waters() = std::vector<Water>();
    signal_modified_hydration_layer();

    // move protein to center of mass
    center();

    // create the grid and hydrate it
    if (grid == nullptr) {create_grid();}
    else {grid->clear_waters();}
    get_waters() = grid->hydrate();
}

ScatteringHistogram Protein::get_histogram() const {
    return phm->calculate_all();
}

Histogram Protein::get_total_histogram() const {
    return phm->calculate();
}

std::shared_ptr<grid::Grid> Protein::get_grid() {
    return grid == nullptr ? create_grid() : grid;
}

void Protein::set_grid(const grid::Grid& grid) {
    this->grid = std::make_shared<grid::Grid>(grid);
}

void Protein::clear_grid() {
    grid = nullptr;
}

void Protein::clear_hydration() {
    if (grid != nullptr) {grid->clear_waters();} // also clear the waters from the grid
    hydration_atoms.clear();
    signal_modified_hydration_layer();
}

unsigned int Protein::body_size() const {
    return bodies.size();
}

unsigned int Protein::atom_size() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0, [] (unsigned int sum, const Body& body) {return sum + body.get_atoms().size();});
}

unsigned int Protein::water_size() const {
    return hydration_atoms.size();
}

Water& Protein::get_waters(unsigned int i) {return hydration_atoms[i];}

const Water& Protein::get_water(unsigned int i) const {return hydration_atoms[i];}

std::vector<double> Protein::calc_debye_scattering_intensity() const {
    std::vector<Atom> atoms = get_atoms();
    const Axis& debye_axis = Axis(settings::axes::bins, settings::axes::qmin, settings::axes::qmax);
    std::vector<double> Q = std::vector<double>(debye_axis.bins);
    double debye_width = debye_axis.width();
    for (unsigned int i = 0; i < debye_axis.bins; i++) {
        Q[i] = debye_axis.min + i*debye_width;
    }

    std::vector<double> I;
    I.reserve(Q.size());
    for (const auto& q : Q) {
        double sum = 0;
        for (const auto& atom_i : atoms) {
            for (const auto& atom_j : atoms) {
                double fi = atom_i.get_effective_charge()*atom_i.get_occupancy();
                double fj = atom_j.get_effective_charge()*atom_j.get_occupancy();
                double qr = q*atom_i.distance(atom_j);
                if (qr < 1e-9) {
                    sum += fi*fj;
                } else {
                    sum += fi*fj*sin(qr)/qr;
                }
            }
        }
        I.push_back(sum);
    }
    return I;
}

void Protein::update_effective_charge(double scaling) {
    static double previous_charge = 0;

    double displaced_vol = scaling*get_volume_grid();
    double displaced_charge = constants::charge::density::water*displaced_vol - previous_charge;
    previous_charge += displaced_charge;

    // number of atoms
    unsigned int N = atom_size();
    double charge_per_atom = -displaced_charge/N;
    if (settings::general::verbose) {
        std::cout << "Total displaced charge: " << displaced_charge << std::endl;
        std::cout << "Added " << charge_per_atom << " electrons to each atom (N: " << N << ")." << std::endl;
    }

    // subtract the charge from all protein atoms
    for (auto& body : bodies) {
        body.update_effective_charge(charge_per_atom);
    }

    updated_charge = true;
}

void Protein::center() {
    if (!centered && settings::protein::center) {
        translate(-get_cm());
        centered = true;
    }
}

void Protein::signal_modified_hydration_layer() const {
    if (phm == nullptr) [[unlikely]] {throw except::unexpected("Protein::signal_modified_hydration_layer: Somehow the histogram manager has not been initialized.");}
    phm->signal_modified_hydration_layer();
}

void Protein::bind_body_signallers() {
    if (phm == nullptr) [[unlikely]] {throw except::unexpected("Protein::bind_body_signallers: Somehow the histogram manager has not been initialized.");}
    for (unsigned int i = 0; i < bodies.size(); i++) {
        bodies[i].register_probe(phm->get_probe(i));
    }
}

std::shared_ptr<fitter::Fit> Protein::fit(const io::ExistingFile& measurement) const {
    hist::ScatteringHistogram h = get_histogram();
    fitter::HydrationFitter fitter(measurement, h);
    return fitter.fit();
}

std::shared_ptr<HistogramManager> Protein::get_histogram_manager() const {return phm;}

void Protein::set_histogram_manager(std::unique_ptr<hist::HistogramManager> manager) {
    phm = std::move(manager);
    bind_body_signallers();
}

void Protein::generate_unit_cell() {
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

void Protein::remove_disconnected_atoms(unsigned int min) {
    if (grid == nullptr) {create_grid();}
    auto to_remove = grid->remove_disconnected_atoms(min);

    // sanity check
    if (to_remove.size() != get_atoms().size()) {
        throw except::unexpected("Protein::remove_disconnected_atoms: "
        "The number of atoms to remove (" + std::to_string(to_remove.size()) + ") does not match the number of protein atoms (" + std::to_string(get_atoms().size()) + ").");
    }

    // remove the atoms from the protein bodies
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

Body& Protein::get_body(unsigned int index) {return bodies[index];}
const Body& Protein::get_body(unsigned int index) const {return bodies[index];}

std::vector<Body>& Protein::get_bodies() {return bodies;}

const std::vector<Body>& Protein::get_bodies() const {return bodies;}