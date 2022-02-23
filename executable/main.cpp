#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <fstream>
#include <iostream>
#include <memory>

using std::string, std::cout, std::endl;

namespace CCP4 {
    /**
     * @brief The 1024-byte header of CCP4 files.
     *        It assumes the endian of the data is the same as the system. 
     */
    struct Header {
        // members CANNOT be reordered!
        int nx;
        int ny;
        int nz;
        int mode;
        int nxstart;
        int nystart;
        int nzstart;
        int mx;
        int my;
        int mz;
        float cella_x;
        float cella_y;
        float cella_z;
        float cellb_alpha;
        float cellb_beta;
        float cellb_gamma;
        int mapc;
        int mapr;
        int maps;
        float dmin;
        float dmax;
        float dmean;
        int ispg;
        int nsymbt;
        std::array<char, 8> extra1;
        std::array<char, 4> exttyp;
        int nversion;
        std::array<char, 84> extra2;
        float origin_x;
        float origin_y;
        float origin_z;
        std::array<char, 4> map;
        unsigned int machst;
        float rms;
        int nlabl;
        std::array<char, 800> label;
    };
    static_assert(sizeof(Header) == 1024, "Size of CCP4 is wrong.");

    
}

using std::vector;

enum class Type{int8, int16, uint16, float16, float32, complex64};

// template<typename T>
// vector<T> create_vector(Type type) {
//     switch (type) {
//         case Type::int8:
//             return vector<short int>();
//         case Type::int16:
//             return vector<int>();
//         case Type::uint16:
//             return vector<unsigned int>();
//         case Type::float16:
//             // return vector<short float>();
//         case Type::float32:
//             return vector<float>();
//         default:
//             //! throw error
//     }
// }

#include <TCanvas.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TStyle.h>
#include "Exceptions.h"

template<typename T>
class Image {
    public:
        Image(const CCP4::Header& header, std::ifstream& istream) : header(header) {
            read(istream);
        }

        void plot(unsigned int layer) const {
            gStyle->SetPalette(kThermometer);
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);

            std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "canvas", 600, 600);
            std::unique_ptr<TH2D> hist = std::make_unique<TH2D>("hist", "hist", header.nx, 0, header.nx*header.cella_x, header.ny, 0, header.ny*header.cella_y);

            double xmin = -3;
            double xmax = 3;
            unsigned int size = 100;
            double step = (xmax - xmin)/size;
            unsigned int err = 0; // number of values which doesn't fit in the histogram
            vector<T> bins(size);
            for (int y = 0; y < header.ny; y++) {
                for (int x = 0; x < header.nx; x++) {
                    unsigned int i = (index(x, y, layer) - xmin)/step;
                    std::cout << "equation: (" << index(x, y, layer) << " - " << xmin << ")/" << size << std::endl;
                    std::cout << i << std::endl;
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

            double negative_limit = 0;
            for (int y = 0; y < header.ny; y++) {
                for (int x = 0; x < header.nx; x++) {
                    T val = index(x, y, layer);
                    bool skip = false;
                    for (const auto range : skip_ranges) {
                        if (range[0] < val && val < range[1]) {
                            hist->SetBinContent(x, y, -1e6);
                            skip = true;
                            break;
                        }
                    }
                    if (skip) {continue;}

                    negative_limit = std::min(negative_limit, double(val));
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

            hist->Draw("cont4z1");

            canvas->SetLeftMargin(0.14);
            canvas->SetRightMargin(0.16);
            canvas->SaveAs("test.pdf");
        }

        void plot3d() const {
            gStyle->SetPalette(kThermometer);
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);

            std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "canvas", 600, 600);
            std::unique_ptr<TH3D> hist = std::make_unique<TH3D>("hist", "hist", header.nx, 0, 1, header.ny, 0, 1, header.nz, 0, 1);
            for (int z = 0; z < header.nz; z++) {
                for (int y = 0; y < header.ny; y++) {
                    for (int x = 0; x < header.nx; x++) {
                        // int my_bin = x + y*header.ny + z*header.nz; // the logical choice of bins
                        // int root_bin = (x + 1) + (y + 1)*(header.ny + 2) + (z + 1); // index 0 is underflow bin and index bins+1 is overflow bin
                        hist->SetBinContent(x, y, z, index(x, y, z));
                    }
                }
            }

            hist->Draw("box");
            canvas->SaveAs("test3d.pdf");
        }

    private:
        const CCP4::Header header;
        // vector<T> data;
        vector<vector<vector<T>>> data;

        void read(std::ifstream& istream) {
            data = vector<vector<vector<T>>>(header.nx, vector<vector<T>>(header.ny, vector<T>(header.nz)));
            for (int i = 0; i < header.nx; i++) {
                for (int j = 0; j < header.ny; j++) {
                    for (int k = 0; k < header.nz; k++) {
                        istream.read(reinterpret_cast<char*>(&data[i][j][k]), sizeof(T));
                    }
                }
            }
        }

        T index(unsigned int x, unsigned int y, unsigned int z) const {
            return data[x][y][z];
        }
};

int main(int argc, char const *argv[]) {
    std::ifstream input("data/A2M_map.ccp4", std::ios::binary);
    CCP4::Header header;
    input.read(reinterpret_cast<char*>(&header), sizeof(header));

    switch(header.mode) {
        case 0: { // int8 --> short int
            Image<short int> image(header, input);
            image.plot(std::stoi(argv[1]));
            break;
        }
        case 1: { // int16 --> int
            Image<int> image(header, input);
            image.plot(std::stoi(argv[1]));
            break;
        }
        case 2: { // float32 --> float
            Image<float> image(header, input);
            image.plot(std::stoi(argv[1]));
            // image.plot3d();
            break;
        }
        case 6: { // uint16 --> unsigned int
            Image<unsigned int> image(header, input);
            image.plot(std::stoi(argv[1]));
            break;
        }
        case 12: { // float16 --> short float
            throw except::parse_error("Short float data format is not currently supported.");
        }
        case 3: // complex32 (2x 16bit int)
            throw except::parse_error("Complex data format is not currently supported.");
        case 4: // complex64 (2x 32bit float)
            throw except::parse_error("Complex data format is not currently supported.");
        default:
            throw except::parse_error("Unrecognized data format (mode " + std::to_string(header.mode) + ").");
    }
    std::cout << header.nx << ", " << header.ny << ", " << header.nz << ", " << header.mode << std::endl;
    std::cout << header.cella_x << ", " << header.cella_y << ", " << header.cella_z << std::endl;
    std::cout << header.nlabl << std::endl;
    std::cout << string(header.label.data()).substr(0, 80) << std::endl;
    return 0;
}