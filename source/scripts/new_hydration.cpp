// includes
#include <vector>
#include <string>

// ROOT
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1.h>

// my own includes
#include "Atom.cpp"
#include "Protein.cpp"
#include "Tools.cpp"

using namespace ROOT;

int main(int argc, char const *argv[])
{
    Protein protein(argv[1]);
    protein.generate_new_hydration();
    protein.save(argv[2]);
    return 0;
}