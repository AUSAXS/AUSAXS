// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_sasview.h>
#include <settings/All.h>
#include <dataset/SimpleDataset.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <hist/detail/SimpleExvModel.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <fitter/SmartFitter.h>
#include <fitter/FitReporter.h>
#include <constants/Constants.h>
#include <utility/Utility.h>

using namespace ausaxs;
using namespace ausaxs::data;

void test_integration(int* test_value) {
    std::cout << "AUSAXS: Starting method \"test_integration\"." << std::endl;
    *test_value += 1;
}

struct {
    std::unique_ptr<data::Molecule> protein;
    std::unique_ptr<SimpleDataset> data;
} iterative_fit_state;

// void iterative_fit_start(
//     double* _data_q, double* _data_I, double* _data_Ierr, int _n_data,
//     double* _pdb_x,  double* _pdb_y,  double* _pdb_z, 
//     const char** _atom_names, const char** _residue_names, const char** _elements, 
//     int _n_pdb, int* _return_status
// ) {
//     std::cout << "AUSAXS: Starting method \"iterative_fit::start\"." << std::endl;

//     // default state is error since we don't trust the input enough to assume success
//     *_return_status = 1;

//     // use the multithreaded version of the simple histogram manager
//     settings::exv::exv_method = settings::exv::ExvMethod::GridScalable;
//     settings::fit::fit_excluded_volume = true;

//     // set qmax as high as it can go
//     settings::axes::qmax = 1;

//     // convert C data
//     std::vector<std::string> atom_names(_n_pdb), residue_names(_n_pdb), elements(_n_pdb);
//     {
//         std::vector<double> q(_data_q, _data_q+_n_data);
//         std::vector<double> I(_data_I, _data_I+_n_data);
//         std::vector<double> Ierr(_data_Ierr, _data_Ierr+_n_data);
//         iterative_fit_state.data = std::make_unique<SimpleDataset>(std::move(q), std::move(I), std::move(Ierr));
//         for (int i = 0; i < _n_pdb; ++i) {
//             atom_names[i]    = std::string(_atom_names[i]);
//             residue_names[i] = std::string(_residue_names[i]);
//             elements[i]      = std::string(_elements[i]);
//         }
//     }

//     std::vector<data::AtomFF> atoms(_n_pdb);
//     for (int i = 0; i < _n_pdb; ++i) {
//         auto atom = constants::symbols::parse_element_string(elements[i]);
//         auto group = constants::symbols::get_atomic_group(residue_names[i], atom_names[i], atom);
//         atoms[i] = data::AtomFF({_pdb_x[i], _pdb_y[i], _pdb_z[i]}, form_factor::get_type(atom, group));
//     }

//     *_return_status = 2;
//     iterative_fit_state.protein = std::make_unique<Molecule>(std::vector{Body{atoms}});
//     iterative_fit_state.protein->generate_new_hydration();
//     *_return_status = 0;
// }

// void iterative_fit_step(double* pars, double* return_I, int* return_status) {
//     std::cout << "AUSAXS: Starting method \"iterative_fit::step\"." << std::endl;
//     std::cout << "DEBUG VERSION!" << std::endl;

//     // default state is error since we don't trust the input enough to assume success
//     *return_status = 1;
//     if (!iterative_fit_state.protein || !iterative_fit_state.data) {return;}

//     *return_status = 2;
//     auto hist = iterative_fit_state.protein->get_histogram();

//     *return_status = 3;
//     double c = pars[0];
//     double d = pars[1];
//     hist->apply_water_scaling_factor(c);
//     static_cast<hist::CompositeDistanceHistogramFFGridSurface*>(hist.get())->apply_excluded_volume_scaling_factor(d);

//     *return_status = 4;
//     auto I = hist->debye_transform(iterative_fit_state.data->x());

//     // write the fitted intensity to the output array
//     *return_status = 5;
//     for (int i = 0; i < static_cast<int>(I.size()); ++i) {
//         return_I[i] = I.y(i);
//     }    
//     *return_status = 0;
// }

// void iterative_fit_finish(double* pars, double* return_I, int* return_status) {
//     std::cout << "AUSAXS: Starting method \"iterative_fit::finish" << std::endl;
//     iterative_fit_step(pars, return_I, return_status);
//     if (*return_status != 0) {return;}
//     *return_status = 1;
//     iterative_fit_state.protein->save("fitted_model.pdb");
//     *return_status = 0;
// }

// void fit_saxs(
//     double* _data_q, double* _data_I, double* _data_Ierr, int _n_data,
//     double* _pdb_x,  double* _pdb_y,  double* _pdb_z, 
//     const char** _atom_names, const char** _residue_names, const char** _elements, 
//     int _n_pdb,
//     double* _return_I, int* _return_status
// ) {
//     std::cout << "AUSAXS: Starting method \"fit_saxs\"." << std::endl;

//     // default state is error since we don't trust the input enough to assume success
//     *_return_status = 1;

//     // use the multithreaded version of the simple histogram manager
//     settings::exv::exv_method = settings::exv::ExvMethod::Simple;
//     settings::fit::fit_excluded_volume = false;

//     // set qmax as high as it can go
//     settings::axes::qmax = 1;

//     // convert C data
//     SimpleDataset data;
//     std::vector<std::string> atom_names(_n_pdb), residue_names(_n_pdb), elements(_n_pdb);
//     {
//         std::vector<double> q(_data_q, _data_q+_n_data);
//         std::vector<double> I(_data_I, _data_I+_n_data);
//         std::vector<double> Ierr(_data_Ierr, _data_Ierr+_n_data);
//         data = SimpleDataset({std::move(q), std::move(I), std::move(Ierr)});
//         for (int i = 0; i < _n_pdb; ++i) {
//             atom_names[i]    = std::string(_atom_names[i]);
//             residue_names[i] = std::string(_residue_names[i]);
//             elements[i]      = std::string(_elements[i]);
//         }
//     }

//     std::vector<data::AtomFF> atoms(_n_pdb);
//     for (int i = 0; i < _n_pdb; ++i) {
//         auto atom = constants::symbols::parse_element_string(elements[i]);
//         auto group = constants::symbols::get_atomic_group(residue_names[i], atom_names[i], atom);
//         atoms[i] = data::AtomFF({_pdb_x[i], _pdb_y[i], _pdb_z[i]}, form_factor::get_type(atom, group));
//     }

//     // construct a molecule from the collection of atom
//     *_return_status = 2;
//     Molecule protein({Body{atoms}});
//     protein.generate_new_hydration();

//     // perform the fit
//     *_return_status = 3;
//     fitter::SmartFitter fitter(std::move(data), protein.get_histogram());
//     auto res = fitter.fit();
//     fitter::FitReporter::report(res.get());
//     fitter::FitReporter::save(res.get(), settings::general::output + "ausaxs_fit_result.txt");

//     res->curves.select_columns({0, 1, 2, 3}).save(
//         settings::general::output + "ausaxs.fit", 
//         "chi2=" + std::to_string(res->fval/res->dof) + " dof=" + std::to_string(res->dof)
//     );

//     // write the fitted intensity to the output array
//     *_return_status = 4;
//     auto fitted_I = res->curves.col(3);
//     for (int i = 0; i < static_cast<int>(fitted_I.size()); ++i) {
//         _return_I[i] = fitted_I[i];
//     }
//     *_return_status = 0;
// }

void debye_no_ff(double* _q, double* _x, double* _y, double* _z, double* _w, int _nq, int _nc, double* _return_Iq, int* _return_status) {
    std::cout << "AUSAXS: Starting method \"evaluate_sans_debye\"." << std::endl;
    // default state is error since we don't trust the input enough to assume success
    *_return_status = 1;

    // use the multithreaded version of the simple histogram manager
    settings::exv::exv_method = settings::exv::ExvMethod::Simple;

    // do not subtract the solvent charge from the atoms
    hist::detail::SimpleExvModel::disable();

    // do not subtract the charge of bound hydrogens
    settings::molecule::implicit_hydrogens = false;

    // set qmax as high as it can go
    settings::axes::qmax = 1;

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