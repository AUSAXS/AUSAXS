#include <em/image.h>

#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TStyle.h>

using namespace em;

Image::Image(std::shared_ptr<ccp4::Header> header) : header(header), data(header->nx, vector<float>(header->ny)) {}

float Image::index(unsigned int x, unsigned int y) const {return data[x][y];}
float& Image::index(unsigned int x, unsigned int y) {return data[x][y];}

vector<Atom> Image::generate_atoms(double cutoff) const {
    vector<Atom> atoms;
    atoms.reserve(128);
    for (int x = 0; x < header->nx; x++) {
        for (int y = 0; y < header->ny; y++) {
            float val = index(x, y);
            if (val < cutoff) {
                continue;
            }

            Vector3 coords{x*header->cella_x, y*header->cella_y, 0};
            atoms.push_back(Atom(0, "C", "", "LYS", "", 0, "", coords, val, 0, "C", ""));
        }
    }

    return atoms;
}

std::unique_ptr<TH2D> Image::as_hist() const {
    std::unique_ptr<TH2D> hist = std::make_unique<TH2D>("hist", "hist", header->nx, 0, header->nx*header->cella_x, header->ny, 0, header->ny*header->cella_y);
    for (int x = 0; x < header->nx; x++) {
        for (int y = 0; y < header->ny; y++) {
            hist->SetBinContent(x, y, index(x, y));
        }
    }
    return hist;
}

void Image::plot_without_solution() const {
    gStyle->SetPalette(kThermometer);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "canvas", 600, 600);
    std::unique_ptr<TH2D> hist = std::make_unique<TH2D>("hist", "hist", header->nx, 0, header->nx*header->cella_x, header->ny, 0, header->ny*header->cella_y);

    float xmin = -3;
    float xmax = 3;
    unsigned int size = 100;
    float step = (xmax - xmin)/size;
    unsigned int err = 0; // number of values which doesn't fit in the histogram
    vector<int> bins(size);
    for (int y = 0; y < header->ny; y++) {
        for (int x = 0; x < header->nx; x++) {
            unsigned int i = (index(x, y) - xmin)/step;
            if (i < 0 || size < i) {
                err++;
            } else {
                bins[i]++;
            }
        }
    }

    // find 2 highest elements in the bin sequence. if they account for more than 10% of the values, we want to remove them visually from the plot.
    int N = header->nx*header->ny;
    vector<vector<float>> skip_ranges;
    unsigned int i1 = std::max_element(bins.begin(), bins.end()) - bins.begin();
    if (bins[i1] > 0.1*N) {skip_ranges.push_back({xmin+i1*step, xmin+(i1+1)*step});}
    bins[i1] = 0;
    unsigned int i2 = std::max_element(bins.begin(), bins.end()) - bins.begin();
    if (bins[i2] > 0.1*N) {skip_ranges.push_back({xmin+i2*step, xmin+(i2+1)*step});}

    // for (unsigned int i = 0; i < bins.size(); i++) {
    //     std::cout << std::to_string(xmin + i*(xmax-xmin)/size) << ": " << bins[i] << std::endl;
    // }

    float negative_limit = 0;
    for (int y = 0; y < header->ny; y++) {
        for (int x = 0; x < header->nx; x++) {
            float val = index(x, y);
            bool skip = false;
            for (const auto range : skip_ranges) {
                if (range[0] < val && val < range[1]) {
                    hist->SetBinContent(x, y, -1e6);
                    skip = true;
                    break;
                }
            }
            if (skip) {continue;}

            negative_limit = std::min(negative_limit, val);
            hist->SetBinContent(x, y, val);
        }
    }
    hist->SetMinimum(negative_limit);

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

    hist->Draw("colz1");

    canvas->SetLeftMargin(0.14);
    canvas->SetRightMargin(0.16);
    canvas->SaveAs("test.pdf");
}