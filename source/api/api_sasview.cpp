// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_sasview.h>

#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/All.h>


using namespace ausaxs;
using namespace ausaxs::data;

void test_integration(int* test_value) {
    *test_value += 1;
}

namespace {
    // the settings debye_no_ff needs, restored to the caller's values when it returns (or throws)
    struct RawDebyeScope {
        settings::exv::ExvMethod exv_method = settings::exv::exv_method;
        bool implicit_hydrogens = settings::molecule::implicit_hydrogens;
        double qmax = settings::axes::qmax;

        RawDebyeScope() {
            settings::exv::exv_method = settings::exv::ExvMethod::None; // simple histogram manager, no solvent charge subtracted from the atoms
            settings::molecule::implicit_hydrogens = false;             // do not add the charge of bound hydrogens
            settings::axes::qmax = 1;                                   // as high as it can go
        }
        ~RawDebyeScope() {
            settings::exv::exv_method = exv_method;
            settings::molecule::implicit_hydrogens = implicit_hydrogens;
            settings::axes::qmax = qmax;
        }
    };
}

void debye_no_ff(double* _q, double* _x, double* _y, double* _z, double* _w, int _nq, int _nc, double* _return_Iq, int* _return_status) {
    // default state is error since we don't trust the input enough to assume success
    *_return_status = 1;
    RawDebyeScope scope;

    // convert coordinate input to Atom objects
    std::vector<double> q(_q, _q+_nq);
    std::vector<data::Atom> atoms(_nc);
    for (int i = 0; i < _nc; ++i) {
        atoms[i] = data::Atom({_x[i], _y[i], _z[i]}, _w[i]);
    }

    // construct a protein from the collection of atom
    *_return_status = 2;
    Molecule protein({Body{atoms}});

    // calculate the distance histogram for the protein
    *_return_status = 3;
    auto dist = protein.get_histogram();

    // perform the Debye transform
    *_return_status = 4;
    auto Iq = dist->debye_transform<false>(q);

    // sanity check - the number of q values should match the number of I(q) values
    if (Iq.size() != _nq) {
        *_return_status = 5;
        return;
    }

    for (int i = 0; i < Iq.size(); ++i) {
        _return_Iq[i] = Iq.y(i);
    }
    *_return_status = 0;
}