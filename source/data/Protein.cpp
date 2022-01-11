#include "data/Protein.h"
#include "data/Body.h"
#include "io/File.h"
#include "Histogram.h"

Protein::Protein(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms) : hydration_atoms(hydration_atoms) {
    bodies = {Body(protein_atoms, this->hydration_atoms)}; // 'this' keyword is necessary, otherwise the objects are bound to the argument instead of the member
    phm = std::make_unique<PartialHistogramManager>(bodies, this->hydration_atoms);
}

Protein::Protein(const vector<vector<Atom>>& protein_atoms, const vector<Hetatom>& hydration_atoms) : hydration_atoms(hydration_atoms) {
    for (size_t i = 0; i < protein_atoms.size(); i++) {
        bodies.push_back(Body(protein_atoms[i], {}));
    }
    phm = std::make_unique<PartialHistogramManager>(bodies, this->hydration_atoms);
}

Protein::Protein(const string& input) {
    bodies = {Body(input)};
    phm = std::make_unique<PartialHistogramManager>(bodies, hydration_atoms);
}

Protein::Protein(const vector<string>& input) {
    for (size_t i = 0; i < input.size(); i++) {
        bodies.push_back(Body(input[i]));
    }
    phm = std::make_unique<PartialHistogramManager>(bodies, hydration_atoms);
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

void Protein::create_grid() {
    grid = std::make_shared<Grid>(setting::grid::base_point, setting::grid::width, setting::grid::bins/setting::grid::width); 
    for (auto const& body : bodies) {
        grid->add(body.protein_atoms);
        // grid->add(body.hydration_atoms);
    }
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
    create_grid();
    hydration_atoms = grid->hydrate();
}

void Protein::calc_histogram() {
    if (!updated_charge && setting::protein::use_effective_charge) {
        update_effective_charge(); // update the effective charge of all proteins. We have to do this since it affects the weights. 
    }
    histogram = std::make_shared<ScatteringHistogram>(phm->calculate_all());
}

shared_ptr<ScatteringHistogram> Protein::get_histogram() {
    if (histogram == nullptr) {calc_histogram();}
    return histogram;
}

void Protein::update_effective_charge() { 
    double displaced_vol = get_volume_grid();
    // double displaced_vol = get_volume_acids();
    double displaced_charge = constants::charge::density::water*displaced_vol;
    // cout << "Displaced volume: " << displaced_vol << ", displaced charge: " << displaced_charge << endl;

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