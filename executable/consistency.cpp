#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;

    // check that we have at least one argument
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mapfile>" << std::endl;
        return 1;
    }

    // load the input map file
    string mapfile = argv[1];
    em::ImageStack map(mapfile); 
    string path = "figures/consistency/" + utility::stem(mapfile) + "/";

    unsigned int evals = 100;
    Dataset data({"dof", "chi2", "cutoff"});
    for (unsigned int i = 0; i < evals; i++) {
        std::cout << "Starting iteration " << i << " of " << evals << std::endl;
        // generate a pdb file from the map at some cutoff
        auto protein = map.get_protein(map.level(3));

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
    }
    data.save(path + "out.txt");
}