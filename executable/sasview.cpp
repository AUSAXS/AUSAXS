#include <settings/All.h>
#include <dataset/Dataset.h>
#include <dataset/SimpleDataset.h>
#include <data/record/Atom.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <plots/PlotHistogram.h>
#include <utility/Utility.h>

#include <fstream>

using namespace data;
using namespace data::record;

// args: <data file> <q file> <output file>
// data file format:   | x | y | z | w |
// q file format:      | q |
// output file format: | q | I(q) |
int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    settings::hist::weighted_bins = true;
    settings::axes::qmax = 1;

    // check args
    if (argc != 4) {
        std::cout << "Wrong number of arguments. See the documentation for help." << std::endl;
        return 1;
    }

    // read files
    std::ifstream data_file(argv[1]);
    std::ifstream q_file(argv[2]);
    if (!data_file.is_open() || !q_file.is_open()) {
        return 1;
    }

    // read data
    std::vector<std::vector<double>> data;
    std::vector<double> q;
    data.reserve(1000);
    q.reserve(1000);
    std::string line;
    while (std::getline(data_file, line)) {
        auto vals = utility::split(line, ' ');
        data.push_back(std::vector<double>{std::stod(vals[0]), std::stod(vals[1]), std::stod(vals[2]), std::stod(vals[3])});
    }
    while (std::getline(q_file, line)) {
        q.push_back(std::stod(line));
    }

    std::vector<Atom> atoms; atoms.reserve(data.size());
    for (unsigned int i = 0; i < data.size(); ++i) {
        atoms.push_back(Atom({data[i][0], data[i][1], data[i][2]}, data[i][3], constants::atom_t::dummy, "", i));
    }

    Molecule protein(atoms);
    auto dist = protein.get_histogram();
    plots::PlotHistogram::quick_plot(*dist.get(), {}, "histogram.png");
    auto Iq = protein.get_histogram()->debye_transform();
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        Iq[i] /= std::exp(-constants::axes::q_vals[i]*constants::axes::q_vals[i]);
    }
    Iq.as_dataset().interpolate(q).save(argv[3]);
    return 1;
}