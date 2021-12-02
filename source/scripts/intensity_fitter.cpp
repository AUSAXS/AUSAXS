// includes
#include <vector>
#include <string>
#include <iostream>

#include "Protein.h"
#include "fitter/IntensityFitter.cpp"

using std::cout, std::endl;

int main(int argc, char const *argv[]) {
    setting::grid::psc = setting::grid::RadialStrategy;

    Protein protein(argv[1]);
    protein.generate_new_hydration();
    std::shared_ptr<Distances> d = protein.get_distances();
    d->set_axes({60, 0, 60});
    auto I = d->calc_debye_scattering_intensity();
    std::vector<double> q = d->get_xaxis();

    IntensityFitter fitter(argv[2], q, I);
    Fitter::Fit result = fitter.fit();
    
    cout << "Result is " << result.params["c"] << "." << endl;
    return 0;
}