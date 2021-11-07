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
#include "plot_style.cpp"

using namespace ROOT;

int main(int argc, char const *argv[])
{
    Protein protein(argv[1]);  
    auto[dp, dh] = protein.calc_distances();
    // protein.generate_new_hydration();
    protein.save("temp2.pdb");
    Protein protein2("temp2.pdb");

    // plots
    setup_style();
    TCanvas* c1 = new TCanvas("c1", "c", 600, 600);
    TH1D* h1 = new TH1D("hist", "h", 500, 0, 50);

    for (int i = 0; i < dp.size(); i++) {
        h1->Fill(dp[i]);
    }
    for (int i = 0; i < dh.size(); i++) {
        h1->Fill(dh[i]);
    }

    h1->Draw("colz");

    std::string path = "figures/main.pdf"; 
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.15);
    c1->SaveAs(path.c_str());


    // save(dp, "temp/distances_protein.txt");
    // save(dh, "temp/distances_hydration.txt");

    // std::cout << "Protein atoms:" << std::endl;
    // for (Atom* a : protein.protein) {
    //     a->print();
    // }

    // std::cout << "Hydration atoms:" << std::endl;
    // for (Atom* a : protein.hydration) {
    //     a->print();
    // }
    return 0;
}