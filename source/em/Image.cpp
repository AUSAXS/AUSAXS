#include <list>

#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <em/Image.h>

using namespace em;

Image::Image(std::shared_ptr<ccp4::Header> header, unsigned int layer) : header(header), data(header->nx, vector<float>(header->ny)), z(layer) {}

void Image::set_z(unsigned int z) {
    this->z = z;
}

float Image::index(unsigned int x, unsigned int y) const {return data[x][y];}
float& Image::index(unsigned int x, unsigned int y) {return data[x][y];}

list<Atom> Image::generate_atoms(double cutoff) const {
    list<Atom> atoms;

    // determine which comparison function to use
    std::function<bool(double)> compare_positive = [&cutoff] (double val) {return val < cutoff;};
    std::function<bool(double)> compare_negative = [&cutoff] (double val) {return val > cutoff;};
    std::function<bool(double)> compare_func = 0 <= cutoff ? compare_positive : compare_negative;

    // loop through all pixels in this image
    double xscale = header->cella_x/header->nx;
    double yscale = header->cella_y/header->ny;
    double zscale = header->cella_z/header->nz;
    for (int x = 0; x < header->nx; x++) {
        for (int y = 0; y < header->ny; y++) {
            float val = index(x, y);
            if (compare_func(val)) {
                continue;
            }

            Vector3 coords{x*xscale, y*yscale, z*zscale};
            atoms.push_back(Atom(0, "C", "", "LYS", "", 0, "", coords, val, 0, "C", ""));
        }
    }

    return atoms;
}

std::unique_ptr<TH2D> Image::as_hist() const {
    std::unique_ptr<TH2D> hist = std::make_unique<TH2D>("hist", "hist", header->nx, 0, header->cella_x, header->ny, 0, header->cella_y);
    for (int x = 0; x < header->nx; x++) {
        for (int y = 0; y < header->ny; y++) {
            hist->SetBinContent(x, y, index(x, y));
        }
    }
    return hist;
}

double Image::mean() const {
    double sum = 0;
    for (int x = 0; x < header->nx; x++) {
        for (int y = 0; y < header->ny; y++) {
            sum += index(x, y);
        }
    }

    return sum/(header->nx * header->ny);
}

Limit Image::limits() const {
    double min = 1e9, max = -1e9;
    for (int x = 0; x < header->nx; x++) {
        for (int y = 0; y < header->ny; y++) {
            double val = index(x, y);
            min = std::min(min, val);
            max = std::max(max, val);
        }
    }

    return Limit(min, max);
}

vector<Limit> Image::minimum_area(double cutoff) const {
}