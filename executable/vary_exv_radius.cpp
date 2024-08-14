#include "fitter/ExcludedVolumeFitter.h"
#include <CLI/CLI.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <grid/Grid.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <fitter/HydrationFitter.h>
#include <dataset/SimpleDataset.h>
#include <settings/All.h>
#include <plots/All.h>

#include <sstream>

class VariableGrid : public grid::Grid {
    public: 
        using grid::Grid::Grid;

		double get_atomic_radius(constants::atom_t) const override {return ra;}
		double get_hydration_radius() const override {return rh;}
        void set_atomic_radius(double ra) {this->ra = ra;}
        void set_hydration_radius(double rh) {this->rh = rh;}
    
    private:
        double ra = 0, rh = 0;
};

int main(int, char const *[]) {
    settings::general::verbose = false;
    settings::general::output = "output/vary_exv_radius/";
    settings::molecule::throw_on_unknown_atom = false;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit;
    settings::hist::fit_excluded_volume = true;
    settings::validate_settings();

    std::vector<std::pair<std::string, std::string>> input_files = {
        {"SASDJG5", "SASDJG5"},
        {"SASDPP4", "consensus/SASDPP4"},
        {"SASDPQ4", "consensus/SASDPQ4"},
        {"SASDPR4", "consensus/SASDPR4"},
        {"SASDPS4", "consensus/SASDPS4"},
        {"SASDPT4", "consensus/SASDPT4"}
    };

    auto fit_single = [] (const std::pair<std::string, std::string>& input, double r) {
        std::string folder = "data/" + input.second + "/";
        data::Molecule protein(folder + input.first + ".pdb");

        std::unique_ptr<grid::Grid> grid = std::make_unique<VariableGrid>(protein.get_bodies());
        static_cast<VariableGrid*>(grid.get())->set_atomic_radius(r);
        protein.set_grid(std::move(grid));
        protein.generate_new_hydration();

        double chi2;
        if (settings::hist::fit_excluded_volume) {
            fitter::ExcludedVolumeFitter fitter({folder + input.first + ".dat"}, protein.get_histogram());
            chi2 = fitter.fit()->fval/fitter.fit()->dof;
        } else {
            fitter::HydrationFitter fitter(folder + input.first + ".dat", protein.get_histogram());
            chi2 = fitter.fit()->fval/fitter.fit()->dof;
        }
        std::cout << input.first << " " << r << " " << chi2 << std::endl;
        return chi2;
    };

    io::Folder(settings::general::output).create();
    std::ofstream file(settings::general::output + "chi2.txt");
    file << "DATA\t";
    std::vector<std::vector<double>> chi2;
    for (double r = 1.5; r <= 3.5; r += 0.1) {
        file << r << "\t";
        std::vector<double> chi2_r;
        for (const auto& input : input_files) {
            chi2_r.push_back(fit_single(input, r));
        }
        chi2.push_back(chi2_r);
    }

    // add current default value
    file << "default" << "\t";
    std::vector<double> chi2_default;
    for (const auto& input : input_files) {
        std::string folder = "data/" + input.second + "/";
        data::Molecule protein(folder + input.first + ".pdb");
        protein.generate_new_hydration();

        if (settings::hist::fit_excluded_volume) {
            fitter::ExcludedVolumeFitter fitter({folder + input.first + ".dat"}, protein.get_histogram());
            chi2_default.push_back(fitter.fit()->fval/fitter.fit()->dof);
        } else {
            fitter::HydrationFitter fitter(folder + input.first + ".dat", protein.get_histogram());
            chi2_default.push_back(fitter.fit()->fval/fitter.fit()->dof);
        }
    }
    chi2.push_back(chi2_default);
    file << std::endl;

    for (unsigned int i = 0; i < input_files.size(); i++) {
        file << input_files[i].first << "\t";
        for (unsigned int j = 0; j < chi2.size(); j++) {
            file << chi2[j][i] << "\t";
        }
        file << std::endl;
    }
}