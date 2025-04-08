#include <CLI/CLI.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <grid/Grid.h>
#include <fitter/SmartFitter.h>
#include <dataset/SimpleDataset.h>
#include <settings/All.h>
#include <plots/All.h>

using namespace ausaxs;

class VariableGrid : public grid::Grid {
    public: 
        using grid::Grid::Grid;

		double get_atomic_radius(form_factor::form_factor_t) const override {return ra;}
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
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface;
    settings::fit::fit_excluded_volume = true;
    settings::axes::qmax = 1;
    settings::validate_settings();

    std::vector<std::pair<std::string, std::string>> input_files = {
        {"1ubq", "1ubq"},
        {"6lyz", "6lyz"},
        {"SASDJG5", "SASDJG5"},
        {"SASDA92", "SASDA92"},
        {"SASDEL8", "SASDEL8"},
        {"SASDFZ3", "SASDFZ3"},
        {"SASDGD2", "SASDGD2"},
        {"SASDAW3", "SASDAW3"},
        {"SASDHP7", "SASDHP7"},
        {"SASDJQ7", "SASDJQ7"},
        {"SASDKG4", "SASDKG4"},
        {"SASDJY3", "SASDJY3"},
        {"SASDMB5", "rigidbody/SASDMB5"},
        {"SASDP39", "rigidbody/SASDP39"},
        {"SASDT85", "SASDT85"},
        {"SASDTT4", "SASDTT4"}
    };

    auto fit_single = [] (const std::pair<std::string, std::string>& input, double r) {
        std::string folder = "data/" + input.second + "/";
        data::Molecule protein(folder + input.first + ".pdb");

        std::unique_ptr<grid::Grid> grid = std::make_unique<VariableGrid>(protein.get_bodies());
        static_cast<VariableGrid*>(grid.get())->set_atomic_radius(r);
        protein.set_grid(std::move(grid));
        protein.generate_new_hydration();

        fitter::SmartFitter fitter(SimpleDataset(folder + input.first + ".dat"), protein.get_histogram());
        auto res = fitter.fit();
        double chi2 = res->fval / res->dof;
        std::cout << input.first << " " << r << " " << chi2 << std::endl;
        return chi2;
    };

    io::Folder(settings::general::output).create();
    std::ofstream file(settings::general::output + "chi2.txt");
    file << "DATA\t";
    std::vector<std::vector<double>> chi2;
    for (double r = 1.5; r <= 3; r += 0.05) {
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

        fitter::SmartFitter fitter(SimpleDataset(folder + input.first + ".dat"), protein.get_histogram());
        auto res = fitter.fit();
        double chi2 = res->fval / res->dof;
        chi2_default.push_back(chi2);
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