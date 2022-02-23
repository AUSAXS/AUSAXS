#include <memory>

#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <em/image.h>
#include <data/Atom.h>
#include <data/Protein.h>
#include <fitter/SimpleIntensityFitter.h>
#include "plots/PlotIntensityFit.h"
#include "plots/PlotIntensityFitResiduals.h"
#include <Exceptions.h>

using namespace em;

ImageStack::ImageStack(string file) {
    std::ifstream input("data/A2M_map.ccp4", std::ios::binary);
    input.read(reinterpret_cast<char*>(&header), sizeof(header));
    read(input, get_byte_size());
}

ImageStack::ImageStack(const ccp4::Header& header, std::ifstream& istream) : header(header) {
    read(istream, get_byte_size());
}

std::unique_ptr<Protein> ImageStack::create_protein(double cutoff) const {
    vector<Atom> atoms;
    atoms.reserve(header.nx);
    for (int x = 0; x < header.nx; x++) {
        for (int y = 0; y < header.ny; y++) {
            for (int z = 0; z < header.nz; z++) {
                float val = index(x, y, z);
                if (val < cutoff) {
                    continue;
                }

                Vector3 coords{x*header.cella_x, y*header.cella_y, z*header.cella_z};
                atoms.push_back(Atom(coords, val, "C", "C", 0));
            }
        }
    }

    return std::make_unique<Protein>(atoms);
}

std::unique_ptr<Grid> ImageStack::create_grid(double cutoff) const {
    vector<Atom> atoms;
    atoms.reserve(header.nx);
    for (int x = 0; x < header.nx; x++) {
        for (int y = 0; y < header.ny; y++) {
            for (int z = 0; z < header.nz; z++) {
                float val = index(x, y, z);
                if (val < cutoff) {
                    continue;
                }

                Vector3 coords{x*header.cella_x, y*header.cella_y, z*header.cella_z};
                atoms.push_back(Atom(coords, val, "C", "C", 0));
            }
        }
    }

    Protein protein(atoms);
    return std::make_unique<Grid>(atoms);
}

ScatteringHistogram ImageStack::calc_scattering_hist() const {
    std::unique_ptr<Protein> protein = create_protein(0.1);
    protein->generate_new_hydration();
    return protein->get_histogram();
}

void ImageStack::fit(string filename) const {
    SimpleIntensityFitter fitter(filename, calc_scattering_hist());
    std::shared_ptr<Fitter::Fit> result = fitter.fit();

    // Fit plot
    PlotIntensityFit plot_f(fitter);
    plot_f.save(output + "intensity_fit." + setting::figures::format);

    // Residual plot
    PlotIntensityFitResiduals plot_r(fitter);
    plot_r.save(output + "residuals." + setting::figures::format);

    result->print();
    cout << "c is: " << result->params["a"]*protein.get_mass()/pow(constants::radius::electron, 2)*constants::unit::mg/pow(constants::unit::cm, 3) << endl;
}

size_t ImageStack::get_byte_size() const {
    switch(header.mode) {
        case 0: // int8 --> short int
            return sizeof(short int);

        case 1: // int16 --> int
            return sizeof(int);

        case 2: // float32 --> float
            return sizeof(float);

        case 6: // uint16 --> unsigned int
            return sizeof(unsigned int);

        case 12: // float16 --> short float
            throw except::parse_error("Short float data format is not currently supported.");

        case 3: // complex32 (2x 16bit int)
            throw except::parse_error("Complex data format is not currently supported.");

        case 4: // complex64 (2x 32bit float)
            throw except::parse_error("Complex data format is not currently supported.");

        default:
            throw except::parse_error("Unrecognized data format (mode " + std::to_string(header.mode) + ").");
    }
}

void ImageStack::plot(unsigned int layer) const {
    gStyle->SetPalette(kThermometer);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "canvas", 600, 600);
    std::unique_ptr<TH2D> hist = std::make_unique<TH2D>("hist", "hist", header.nx, 0, header.nx*header.cella_x, header.ny, 0, header.ny*header.cella_y);

    for (int x = 0; x < header.nx; x++) {
        for (int y = 0; y < header.ny; y++) {
            hist->SetBinContent(x, y, index(x, y, layer));
        }
    }

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

    hist->Draw("cont4z");

    canvas->SetLeftMargin(0.14);
    canvas->SetRightMargin(0.16);
    canvas->SaveAs("test.pdf");
}

void ImageStack::plot_no_solution(unsigned int layer) const {
    gStyle->SetPalette(kThermometer);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "canvas", 600, 600);
    std::unique_ptr<TH2D> hist = std::make_unique<TH2D>("hist", "hist", header.nx, 0, header.nx*header.cella_x, header.ny, 0, header.ny*header.cella_y);

    float xmin = -3;
    float xmax = 3;
    unsigned int size = 100;
    float step = (xmax - xmin)/size;
    unsigned int err = 0; // number of values which doesn't fit in the histogram
    vector<int> bins(size);
    for (int y = 0; y < header.ny; y++) {
        for (int x = 0; x < header.nx; x++) {
            unsigned int i = (index(x, y, layer) - xmin)/step;
            if (i < 0 || size < i) {
                err++;
            } else {
                bins[i]++;
            }
        }
    }

    // find 2 highest elements in the bin sequence. if they account for more than 10% of the values, we want to remove them visually from the plot.
    int N = header.nx*header.ny;
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
    for (int y = 0; y < header.ny; y++) {
        for (int x = 0; x < header.nx; x++) {
            float val = index(x, y, layer);
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

void ImageStack::plot3d() const {
    gStyle->SetPalette(kThermometer);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "canvas", 600, 600);
    std::unique_ptr<TH3D> hist = std::make_unique<TH3D>("hist", "hist", header.nx, 0, 1, header.ny, 0, 1, header.nz, 0, 1);
    for (int z = 0; z < header.nz; z++) {
        for (int y = 0; y < header.ny; y++) {
            for (int x = 0; x < header.nx; x++) {
                hist->SetBinContent(x, y, z, index(x, y, z));
            }
        }
    }

    hist->Draw("box");
    canvas->SaveAs("test3d.pdf");
}

void ImageStack::read(std::ifstream& istream, size_t byte_size) {
    data = vector<vector<vector<float>>>(header.nx, vector<vector<float>>(header.ny, vector<float>(header.nz)));
    for (int i = 0; i < header.nx; i++) {
        for (int j = 0; j < header.ny; j++) {
            for (int k = 0; k < header.nz; k++) {
                istream.read(reinterpret_cast<char*>(&index(i, j, k)), byte_size);
                // istream.read(reinterpret_cast<char*>(&data[k][i][j]), byte_size);
            }
        }
    }
}

float& ImageStack::index(unsigned int x, unsigned int y, unsigned int z) {
    return data[x][y][z];
}

float ImageStack::index(unsigned int x, unsigned int y, unsigned int z) const {
    return data[x][y][z];
}