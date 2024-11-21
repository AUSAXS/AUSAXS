#include <CLI/CLI.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <fitter/SmartFitter.h>
#include <dataset/SimpleDataset.h>
#include <utility/Utility.h>
#include <settings/All.h>
#include <plots/All.h>

#include <fstream>
#include <iomanip>

using namespace ausaxs;

int main(int, char const*[]) {
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid;
    settings::molecule::use_effective_charge = false;
    settings::fit::fit_exv_debye_waller = true;

    SimpleDataset saxs("data/hub_data/193l_lyz/193l_lyz.dat");
    data::Molecule protein("data/hub_data/193l_lyz/193l_lyz.pdb");
    protein.add_implicit_hydrogens();

    unsigned int N = 1000; 
    std::vector<double> exv_vals, atomic_vals, exv_single_vals;
    for (unsigned int i = 0; i < N; ++i) {
        std::cout << "Progress: " << double(i)/N*100 << "%\r" << std::flush;
        protein.generate_new_hydration();

        settings::fit::fit_atomic_debye_waller = true;
        fitter::SmartFitter fitter(saxs, protein.get_histogram());
        auto result = fitter.fit();
        exv_vals.push_back(   result->get_parameter(constants::fit::Parameters::DEBYE_WALLER_EXV));
        atomic_vals.push_back(result->get_parameter(constants::fit::Parameters::DEBYE_WALLER_ATOMIC));

        settings::fit::fit_atomic_debye_waller = false;
        fitter::SmartFitter fitter_single(saxs, protein.get_histogram());
        auto result_single = fitter_single.fit();
        exv_single_vals.push_back(result_single->get_parameter(constants::fit::Parameters::DEBYE_WALLER_EXV));
    }

    std::ofstream out("temp/debye_waller_fitting.txt");
    out << "exv\tatomic\texv_single\n";
    for (unsigned int i = 0; i < N; ++i) {
        out << std::setprecision(3) << exv_vals[i] << "\t" << std::setprecision(4) << atomic_vals[i] << "\t" << std::setprecision(4) << exv_single_vals[i] << "\n";
    }

    return 0;
}