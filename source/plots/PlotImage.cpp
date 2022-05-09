#include <plots/PlotImage.h>
#include <em/Image.h>
#include <utility/Utility.h>

#include <TH2D.h>
#include <TGraph.h>
#include <THStack.h>
#include <TCanvas.h>

plots::PlotImage::PlotImage(const em::Image& image) : image(image) {
    canvas = std::make_unique<TCanvas>(utility::uid("canvas").c_str(), "canvas", 1200, 1200);
    pad1 = std::make_unique<TPad>("PlotImagePad1", "pad1", 0, 0, 1, 1);
    pad2 = std::make_unique<TPad>("PlotImagePad2", "pad2", 0, 0, 1, 1);

    pad1->SetRightMargin(0.15);
    pad2->SetRightMargin(0.15);
    pad1->SetLeftMargin(0.11);
    pad2->SetLeftMargin(0.11);

    pad1->Draw();
    pad1->cd();
    std::unique_ptr<TH2D> hist = plot_hist();

    pad2->SetFillStyle(0);
    pad2->SetFrameFillStyle(0);
    pad2->Draw();

    // set correct axis ranges on pad2
    Double_t bm = pad1->GetBottomMargin();
    Double_t lm = pad1->GetLeftMargin();
    Double_t rm = pad1->GetRightMargin();
    Double_t to = pad1->GetTopMargin();
    Double_t x1 = hist->GetXaxis()->GetXmin();
    Double_t yf = hist->GetYaxis()->GetXmin();
    Double_t x2 = hist->GetXaxis()->GetXmax();
    Double_t y2 = hist->GetYaxis()->GetXmax();

    Double_t Xa = (x2-x1)/(1-lm-rm)-(x2-x1);
    Double_t Ya = (y2-yf)/(1-bm-to)-(y2-yf);
    Double_t LM = Xa*(lm/(lm+rm));
    Double_t RM = Xa*(rm/(lm+rm));
    Double_t BM = Ya*(bm/(bm+to));
    Double_t TM = Ya*(to/(bm+to));
    pad2->Range(x1-LM,yf-BM,x2+RM,y2+TM);
    pad2->cd();
}

plots::PlotImage::~PlotImage() = default;

void plots::PlotImage::save(std::string path) const {
    utility::create_directories(path);
    canvas->SaveAs(path.c_str());
}

void plots::PlotImage::plot_atoms(double cutoff) const {
    const std::list<Atom>& atoms = image.generate_atoms(cutoff);
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(atoms.size());
    y.reserve(atoms.size());
    for (const Atom& atom : atoms) {
        x.push_back(atom.coords.x());
        y.push_back(atom.coords.y());
    }
    std::unique_ptr<TGraph> graph = std::make_unique<TGraph>(x.size(), x.data(), y.data());
    graph->SetMarkerStyle(kFullDotSmall);
    graph->SetMarkerSize(1.2);
    graph->DrawClone("p");
}

std::unique_ptr<TH2D> plots::PlotImage::plot_hist() const {
    gStyle->SetPalette(kThermometer);

    std::unique_ptr<TH2D> hist = image.as_hist();

    hist->GetXaxis()->SetTitle("Length [Angstrom]");
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetNdivisions(204);

    hist->GetYaxis()->SetTitle("Length [Angstrom]");
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetNdivisions(204);

    hist->GetZaxis()->SetTitle("Electron density [?]");
    hist->GetZaxis()->CenterTitle();
    hist->GetZaxis()->SetTitleOffset(1.3);
    hist->GetZaxis()->SetNdivisions(505);

    hist->DrawClone("cont4z");
    return hist;
}