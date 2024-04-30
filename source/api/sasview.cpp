/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

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
    std::cout << "AUSAXS: Started evaluating Debye equation." << std::endl;
    // default state is error since we don't trust the input enough to assume success
    *_return_status = 1;

    // use the multithreaded version of the simple histogram manager
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;

    // do not subtract the solvent charge from the atoms
    settings::molecule::use_effective_charge = false;

    // do not subtract the charge of bound hydrogens
    settings::molecule::implicit_hydrogens = false;

    // use weighted bins for the histogram approach - this dramatically improves the accuracy
    settings::hist::weighted_bins = true;

    // set qmax as high as it can go. Values beyond this are suppoted, but will recalculate the sinc(x) lookup table at runtime
    settings::axes::qmax = 1;

    // convert C input to C++
    std::vector<double> q(_q, _q+_nq);
    std::vector<double> x(_x, _x+_nc);
    std::vector<double> y(_y, _y+_nc);
    std::vector<double> z(_z, _z+_nc);
    std::vector<double> w(_w, _w+_nc);

    // convert coordinate input to the Atom object
    std::vector<Atom> atoms(_nc);
    for (int i = 0; i < _nc; ++i) {
        atoms[i] = Atom({x[i], y[i], z[i]}, w[i], constants::atom_t::dummy, "", i);
    }

    // construct a protein from the collection of atom
    *_return_status = 2;
    Molecule protein(atoms);

    // calculate the distance histogram for the protein
    *_return_status = 3;
    auto dist = protein.get_histogram();

    // perform the Debye transform
    *_return_status = 4;
    auto Iq = dist->debye_transform(q);

    // sanity check - the number of q values should match the number of I(q) values
    if ((int) Iq.size() != _nq) {
        *_return_status = 5;
        return;
    }

    // remove the form factor applied by the debye transform
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        _return_Iq[i] =  Iq.y(i) / std::exp(-std::pow(Iq.x(i), 2));
    }
    *_return_status = 0;
}