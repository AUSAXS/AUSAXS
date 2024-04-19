/*
This software is distributed under the GNU General Public License v3.0. 
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
    std::cout << "CHECKPOINT 1" << std::endl;

    // default state is error since we don't trust the input enough to assume success
    *_return_status = 1;
    std::cout << "CHECKPOINT 2" << std::endl;

    // use the multithreaded version of the simple histogram manager
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;
    std::cout << "CHECKPOINT 3" << std::endl;

    // do not subtract the solvent charge from the atoms
    settings::molecule::use_effective_charge = false;
    std::cout << "CHECKPOINT 4" << std::endl;

    // do not subtract the charge of bound hydrogens
    settings::molecule::implicit_hydrogens = false;
    std::cout << "CHECKPOINT 5" << std::endl;

    // use weighted bins for the histogram approach - this dramatically improves the accuracy
    settings::hist::weighted_bins = true;
    std::cout << "CHECKPOINT 6" << std::endl;

    // set qmax as high as it can go. Values beyond this are suppoted, but will recalculate the sinc(x) lookup table at runtime
    settings::axes::qmax = 1;
    std::cout << "CHECKPOINT 7" << std::endl;

    // convert C input to C++
    std::cout << "CHECKPOINT 8" << std::endl;
    std::vector<double> q(_q, _q+_nq);
    std::cout << "CHECKPOINT 9" << std::endl;
    std::vector<double> x(_x, _x+_nc);
    std::cout << "CHECKPOINT 10" << std::endl;
    std::vector<double> y(_y, _y+_nc);
    std::cout << "CHECKPOINT 11" << std::endl;
    std::vector<double> z(_z, _z+_nc);
    std::cout << "CHECKPOINT 12" << std::endl;
    std::vector<double> w(_w, _w+_nc);
    std::cout << "CHECKPOINT 13" << std::endl;

    // convert coordinate input to the Atom object
    std::cout << "CHECKPOINT 14" << std::endl;
    std::vector<Atom> atoms(_nc);
    std::cout << "CHECKPOINT 15" << std::endl;
    for (int i = 0; i < _nc; ++i) {
        atoms[i] = Atom({x[i], y[i], z[i]}, w[i], constants::atom_t::dummy, "", i);
    }
    std::cout << "CHECKPOINT 16" << std::endl;

    // construct a protein from the collection of atom
    Molecule protein(atoms);
    std::cout << "CHECKPOINT 17" << std::endl;

    // calculate the distance histogram for the protein
    auto dist = protein.get_histogram();
    std::cout << "CHECKPOINT 18" << std::endl;

    // perform the Debye transform
    auto Iq = dist->debye_transform(q);
    std::cout << "CHECKPOINT 19" << std::endl;

    // sanity check - the number of q values should match the number of I(q) values
    if ((int) Iq.size() != _nq) {
        *_return_status = 2;
        std::cout << "CHECKPOINT 20a" << std::endl;
        return;
    }

    // remove the form factor applied by the debye transform
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        _return_Iq[i] =  Iq.y(i) / std::exp(-std::pow(Iq.x(i), 2));
    }
    std::cout << "CHECKPOINT 20b" << std::endl;
    *_return_status = 0;
}