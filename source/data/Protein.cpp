#include <cassert>

#include <data/Protein.h>
#include <data/Body.h>
#include <io/File.h>
#include <Histogram.h>

Protein::Protein(const vector<Body>& bodies, const vector<Hetatom>& hydration_atoms) : hydration_atoms(hydration_atoms), bodies(bodies) {
    phm = std::make_unique<PartialHistogramManager>(this);
    bind_body_signallers();
}

Protein::Protein(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms) : hydration_atoms(hydration_atoms) {
    bodies = {Body(protein_atoms, this->hydration_atoms)}; // 'this' keyword is necessary, otherwise the objects are bound to the argument instead of the member
    phm = std::make_unique<PartialHistogramManager>(this);
    bind_body_signallers();
}

Protein::Protein(const vector<vector<Atom>>& protein_atoms, const vector<Hetatom>& hydration_atoms) : hydration_atoms(hydration_atoms) {
    for (size_t i = 0; i < protein_atoms.size(); i++) {
        bodies.push_back(Body(protein_atoms[i], vector<Hetatom>(0)));
    }
    phm = std::make_unique<PartialHistogramManager>(this);
    bind_body_signallers();
}

Protein::Protein(Protein&& protein) noexcept : hydration_atoms(std::move(protein.hydration_atoms)), bodies(std::move(protein.bodies)), updated_charge(protein.updated_charge), centered(protein.centered) {
    phm = std::make_unique<PartialHistogramManager>(this);
    bind_body_signallers();
}

Protein::Protein(string input) {
    Body b1(input);
    bodies = {b1};
    hydration_atoms = std::move(bodies[0].hydration_atoms);
    bodies[0].hydration_atoms.clear();
    phm = std::make_unique<PartialHistogramManager>(this);
    bind_body_signallers();
}

Protein::Protein(const vector<string>& input) {
    for (size_t i = 0; i < input.size(); i++) {
        bodies.push_back(Body(input[i]));
    }
    phm = std::make_unique<PartialHistogramManager>(this);
    bind_body_signallers();
}

void Protein::translate(const Vector3& v) {
    for (auto& body : bodies) {
        body.translate(v);
    }
}

void Protein::save(string path) {
    // if there's only a single body, just save that instead
    if (bodies.size() == 1) {
        bodies[0].hydration_atoms = hydration_atoms;
        bodies[0].save(path);
        return;
    }

    // otherwise we'll have to create a new file
    File file(get_protein_atoms(), hydration_atoms);
    file.write(path);
}

double Protein::get_volume_acids() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.get_volume_acids();});
}

double Protein::get_mass() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0.0, [] (double sum, const Body& body) {return sum + body.get_mass();});
}

double Protein::get_volume_grid() {
    if (grid == nullptr) {create_grid();}
    return grid->get_volume();
}

shared_ptr<Grid> Protein::create_grid() {
    grid = std::make_shared<Grid>(bodies); 
    return grid;
}

vector<Atom> Protein::get_protein_atoms() const {
    int N = std::accumulate(bodies.begin(), bodies.end(), 0, [] (double sum, const Body& body) {return sum + body.protein_atoms.size();});
    vector<Atom> atoms(N);
    int n = 0; // current index
    for (const auto& body : bodies) {
        for (const auto& a : body.protein_atoms) {
            atoms[n] = a;
            n++;
        }
    }
    assert(n == N);
    return atoms;
}

Vector3 Protein::get_cm() const {
    Vector3 cm;
    double M = 0; // total mass

    // iterate through all constituent bodies
    for (const auto& body : bodies) {
        // iterate through their protein atoms
        for (const auto& a : body.protein_atoms) {
            double m = a.get_mass();
            M += m;
            cm += a.coords*m;
        }

        // iterate through their hydration atoms
        for (const auto& a : body.hydration_atoms) {
            double m = a.get_mass();
            M += m;
            cm += a.coords*m;
        }
    }

    // iterate through any generated hydration atoms
    for (const auto& a : hydration_atoms) {
        double m = a.get_mass();
        M += m;
        cm += a.coords*m;
    }

    return cm/M;
}

vector<Hetatom> Protein::get_hydration_atoms() const {return hydration_atoms;}

void Protein::generate_new_hydration() {
    // delete the old hydration layer
    hydration_atoms = vector<Hetatom>();
    phm->signal_modified_hydration_layer();

    // move protein to center of mass
    center();

    // create the grid and hydrate it
    if (grid == nullptr) {create_grid();}
    else {grid->clear_waters();}
    hydration_atoms = grid->hydrate();
}

ScatteringHistogram Protein::get_histogram() {
    if (!updated_charge && setting::protein::use_effective_charge) {
        update_effective_charge(); // update the effective charge of all proteins. We have to do this since it affects the weights. 
    }
    return phm->calculate_all();
}

Histogram Protein::get_total_histogram() const {
    return phm->calculate();
}

shared_ptr<Grid> Protein::get_grid() {
    return grid == nullptr ? create_grid() : grid;
}

void Protein::set_grid(const Grid& grid) {
    this->grid = std::make_shared<Grid>(grid);
}

void Protein::clear_grid() {
    grid = nullptr;
}

void Protein::clear_hydration() {
    hydration_atoms.clear();
}

size_t Protein::body_size() const {
    return bodies.size();
}

size_t Protein::atom_size() const {
    return std::accumulate(bodies.begin(), bodies.end(), 0, [] (size_t sum, const Body& body) {return sum + body.protein_atoms.size();});
}

vector<double> Protein::calc_debye_scattering_intensity() {
    if (!updated_charge && setting::protein::use_effective_charge) {
        update_effective_charge(); // update the effective charge of all proteins. We have to do this since it affects the weights. 
    }

    vector<Atom> atoms = get_protein_atoms();
    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;
    vector<double> Q = vector<double>(debye_axis.bins);
    double debye_width = debye_axis.width();
    for (unsigned int i = 0; i < debye_axis.bins; i++) {
        Q[i] = debye_axis.min + i*debye_width;
    }

    vector<double> I;
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

void Protein::update_effective_charge() { 
    double displaced_vol = get_volume_grid();
    // double displaced_vol = get_volume_acids();
    double displaced_charge = constants::charge::density::water*displaced_vol;
    // cout << "Volume: acid: " << get_volume_acids() << ", grid: " << displaced_vol << endl;
    cout << "Displaced charge: " << displaced_charge << endl;

    // number of atoms
    int N = std::accumulate(bodies.begin(), bodies.end(), 0, [] (double sum, const Body& body) {return sum + body.protein_atoms.size();});
    double charge_per_atom = -displaced_charge/N;
    cout << "Added " << charge_per_atom << " additional charge to each protein atom (N: " << N << ")." << endl;

    // subtract the charge from all protein atoms
    for (auto& body : bodies) {
        body.update_effective_charge(charge_per_atom);
    }

    updated_charge = true;
}

void Protein::center() {
    if (!centered && setting::protein::center) {
        translate(-get_cm());
        centered = true;
    }
}

void Protein::bind_body_signallers() {
    if (phm == nullptr) {throw except::unexpected("Error in Protein::bind_body_signallers: Somehow the histogram manager has not been initialized.");}
    for (unsigned int i = 0; i < bodies.size(); i++) {
        std::cout << "Registering probe to body " << i << std::endl;
        bodies[i].register_probe(phm->get_probe(i));
    }
}

std::shared_ptr<PartialHistogramManager> Protein::get_histogram_manager() const {return phm;}