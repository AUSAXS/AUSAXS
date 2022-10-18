#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

#include <filesystem>
#include <fstream>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::plot::em::plot_cutoff_points = true;
    setting::em::sample_frequency = 2;
    setting::em::save_pdb = false;
    setting::em::hydrate = false;

    // check that we have at least one argument
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mapfile>" << std::endl;
        return 1;
    }

    // load the input map file
    string mapfile = argv[1];
    em::ImageStack map(mapfile); 
    string path = "figures/consistency/" + utility::stem(mapfile) + "/";

    unsigned int evals = 1000;
    Dataset data({"dof", "chi2", "cutoff"});
    for (unsigned int i = 0; i < evals; i++) {
        std::cout << "Starting iteration " << i+1 << " of " << evals << std::endl;

        // load all files from previous runs
        if (std::filesystem::exists(path + "fits/" + std::to_string(i) + ".txt")) {
            std::cout << "\tFound existing fit" << std::endl;
            std::fstream file(path + "fits/" + std::to_string(i) + ".txt");
            std::string line;

            // skip the headers
            for (unsigned int j = 0; j < 5; j++) {
                std::getline(file, line);
            }
            auto tokens = utility::split(line, "| ");
            double chi2 = std::stod(tokens[1]);
            double dof = std::stod(tokens[3]);

            // skip stuff until we get to the parameters
            for (unsigned int j = 0; j < 3; j++) {
                std::getline(file, line);
            }

            // read the cutoff
            double cutoff = -1;
            for (unsigned int j = 0; j < 3; j++) {
                tokens = utility::split(line, "| ");            
                if (tokens[0] == "cutoff") {
                    cutoff = std::stod(tokens[1]);
                    break;
                }
            }

            if (cutoff == -1) {
                throw except::io_error("main: Could not find cutoff");
            }
            data.push_back({dof, chi2, cutoff});
            std::cout << "\tSuccessfully read existing fit" << std::endl;
            continue;
        }

        // generate a pdb file from the map at some cutoff
        auto protein = map.get_protein(map.level(3));

        if (setting::em::hydrate) {
            protein->generate_new_hydration();
        }

        // generate a measurement from the protein
        string mfile = path + "test.RSR";
        auto m = protein->get_histogram().calc_debye_scattering_intensity();
        m.reduce(100);
        m.simulate_errors();
        m.simulate_noise();
        m.save(mfile);

        // fit the measurement to the map
        auto res = map.fit(mfile);

        data.push_back({double(res->dof), res->fval, res->get_parameter("cutoff")});
        FitReporter::report(res);
        FitReporter::save(path + "fits/" + std::to_string(i) + ".txt", res);
        plots::PlotIntensityFit::quick_plot(res, path + "intensity_fit.pdf");
        plots::PlotIntensityFitResiduals::quick_plot(res, path + "residuals.pdf");

        // delete the temporary files
        remove(mfile.c_str());
    }
    data.save(path + "out.txt", "true_cutoff: " + std::to_string(map.level(3)));
}