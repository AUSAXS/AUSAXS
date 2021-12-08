// includes
#include <vector>
#include <string>

// ROOT
// #include <TStyle.h>
// #include <TROOT.h>
// #include <TCanvas.h>
// #include <TH1.h>

// my own includes
#include "data/Atom.h"
#include "Protein.h"
#include "Tools.h"
#include "plot_style.h"

using namespace ROOT;

int main(int argc, char const *argv[]) {
    setting::grid::psc = setting::grid::RadialStrategy;
    setting::grid::ra = 3;
    setting::grid::rh = 3;

    Protein protein(argv[1]);
    protein.generate_new_hydration();
    protein.generate_volume_file(argv[2]);
    return 0;
}