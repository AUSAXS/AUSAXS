#include <api/sasview.h>

#include <settings/All.h>
#include <dataset/Dataset.h>
#include <dataset/SimpleDataset.h>
#include <data/record/Atom.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <utility/Utility.h>
#include <utility/MultiThreading.h>

using namespace data;
using namespace data::record;

void evaluate_sans_debye(double* _q, double* _x, double* _y, double* _z, double* _w, int _nq, int _nc, int* _return_status, double* _return_Iq) {
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManager;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    settings::hist::weighted_bins = false;
    settings::axes::qmax = 1;

    std::vector<double> q(_q, _q+_nq);
    std::vector<double> x(_x, _x+_nc);
    std::vector<double> y(_y, _y+_nc);
    std::vector<double> z(_z, _z+_nc);
    std::vector<double> w(_w, _w+_nc);

    *_return_status = 1;
    if (q.front() < constants::axes::q_axis.min || constants::axes::q_axis.max < q.back()) {
        *_return_status = 2;
        return;
    }

    std::vector<Atom> atoms(_nc);
    for (unsigned int i = 0; i < _nc; ++i) {
        atoms[i] = Atom({x[i], y[i], z[i]}, w[i], constants::atom_t::dummy, "", i);
    }

    Molecule protein(atoms);
    auto dist = protein.get_histogram();
    auto Iq = protein.get_histogram()->debye_transform();

    // remove the form factor applied by the debye transform
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        Iq[i] /= std::exp(-constants::axes::q_vals[i]*constants::axes::q_vals[i]);
    }

    auto res = Iq.as_dataset().interpolate(q);
    for (unsigned int i = 0; i < _nq; ++i) {
        _return_Iq[i] = res.y(i);
    }
    *_return_status = 0;

    std::cout << "Successfully calculated I(q) for " << _nc << " atoms and " << _nq << " q values." << std::endl;
}