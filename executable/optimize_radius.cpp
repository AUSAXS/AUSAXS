// #include <vector>
// #include <string>
// #include <float.h>

// #include "data/Protein.h"
// #include "utility/Settings.h"

// #include "TCanvas.h"
// #include "TGraph.h"
// #include "TLine.h"

// class debug_protein : public Body {
// public:
//     using Body::Body;
//     void refresh_grid() {create_grid();}
// };

// int main(int, char const *argv[]) {
//     debug_protein protein(argv[1]);
//     double acid_volume = protein.get_volume_acids();
//     setting::grid::width = 0.3;

//     std::cout << "Average radius is " << cbrt(acid_volume/protein.protein_atoms.size()) << std::endl;
//     std::cout << "acid_vol: " << acid_volume << ", N: " << protein.protein_atoms.size() << std::endl;

//     std::vector<double> ra(20);
//     std::vector<double> diff(ra.size());
//     double mindiff = DBL_MAX;
//     for (size_t i = 0; i < ra.size(); i++) {
//         ra[i] = 1+i*0.1;
//         setting::grid::ra = ra[i];
//         protein.refresh_grid();
//         double vol = protein.get_volume_grid();
//         diff[i] = abs(acid_volume - vol);
//         if (diff[i] < mindiff) {mindiff = diff[i];}
//     }

//     std::unique_ptr<TCanvas> c1 = std::make_unique<TCanvas>("c1", "canvas", 600, 600);
//     std::unique_ptr<TGraph> g = std::make_unique<TGraph>(ra.size(), &ra[0], &diff[0]);
//     g->Draw("AL");

//     // setup the canvas and save the plot
//     string path = string(argv[2]) + "radius_optimizer." + setting::figures::format;
//     c1->SetRightMargin(0.15);
//     c1->SetLeftMargin(0.15);
//     c1->SaveAs(path.c_str());

//     return 0;
// }