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
    Protein protein(argv[1]);
    protein.generate_new_hydration();
    protein.generate_volume_file(argv[2]);
    return 0;
}