#include <crystal/CrystalScattering.h>
#include <crystal/Fval.h>
#include <utility/Exceptions.h>
#include <crystal/miller/AllMillers.h>
#include <crystal/miller/FibonacciMillers.h>
#include <crystal/miller/ReducedMillers.h>
#include <crystal/io/CrystalReaderFactory.h>
#include <crystal/miller/MillerGenerationFactory.h>
#include <crystal/miller/MillerGenerationStrategy.h>
#include <crystal/io/CrystalReader.h>
#include <settings/CrystalSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <hydrate/Grid.h>
#include <dataset/SimpleDataset.h>
#include <utility/Basis3D.h>

#include <atomic>
#include <thread>

using namespace crystal;

CrystalScattering::CrystalScattering(const std::string& input) {
    initialize();
    auto reader = factory::CrystalReaderFactory::create(input);
    auto [bases, points] = reader->read(input);
    // std::cout << bases.x << std::endl;
    // std::cout << bases.y << std::endl;
    // std::cout << bases.z << std::endl;
    Fval::set_points(std::move(points));
    Fval::set_basis(bases);
}

CrystalScattering::CrystalScattering(const grid::Grid& grid) {
    initialize();
    convert_grid(grid);
}

void CrystalScattering::initialize() {
    miller_strategy = factory::construct_miller_strategy();
}

SimpleDataset CrystalScattering::calculate() const {
    if (Fval::get_points().empty()) {throw except::invalid_argument("CrystalScattering::calculate: No points were set.");}
    if (Fval::get_basis().x.x() == 0 || Fval::get_basis().y.y() == 0 || Fval::get_basis().z.z() == 0) {throw except::invalid_argument("CrystalScattering::calculate: No basis was set.");}
    auto millers = miller_strategy->generate();

    std::vector<Fval> fvals(millers.size());
    std::atomic<unsigned int> index = 0;
    auto dispatcher = [&] () {
        while (true) {
            unsigned int start = index.fetch_add(1000);
            unsigned int end = std::min<unsigned int>(start + 1000, millers.size());
            if (start >= millers.size()) {break;}
            std::cout << start << "/" << millers.size() << "          \r" << std::flush;
            for (unsigned int i = start; i < end; i++) {
                fvals[i] = Fval(millers[i].h, millers[i].k, millers[i].l);
            }
        }
    };

    // start threads
    std::vector<std::thread> threads;
    for (unsigned int i = 0; i < settings::general::threads; i++) {
        threads.push_back(std::thread(dispatcher));
    }

    // wait for threads to finish
    for (auto& thread : threads) {
        thread.join();
    }

    // sort fvals by q
    std::sort(fvals.begin(), fvals.end(), [] (const Fval& a, const Fval& b) {return a.qlength < b.qlength;});

    // prepare 100 equidistant logarithmic bins between 0 and 0.5
    std::vector<double> bins(settings::axes::bins);
    double logmin = std::log10(settings::axes::qmin);
    double logmax = std::log10(settings::axes::qmax);
    double logstep = (logmax - logmin)/bins.size();
    for (unsigned int i = 0; i < bins.size(); i++) {
        bins[i] = std::pow(10, logmin + i*logstep);
    }

    // bin the data
    SimpleDataset data;
    unsigned int bin_index = 0;
    for (unsigned int i = 0; i < bins.size()-1; i++) {
        double qmin = bins[i];
        double qmax = bins[i+1];
        double Isum = 0;
        unsigned int count = 0;
        while (bin_index < fvals.size() && fvals[bin_index].qlength < qmax) {
            Isum += fvals[bin_index].I();
            count++;
            bin_index++;
        }
        data.push_back(qmin, count == 0 ? 0 : Isum/count);
    }

    return data;
}

void CrystalScattering::convert_grid(const grid::Grid& grid) const {
    auto axes = grid.get_axes();
    std::vector<Vector3<double>> points(axes.x.bins*axes.y.bins*axes.z.bins);
    unsigned int index = 0;
    double xstep = axes.x.step(), ystep = axes.y.step(), zstep = axes.z.step();
    for (unsigned int i = 0; i < axes.x.bins; i++) {
        for (unsigned int j = 0; j < axes.y.bins; j++) {
            for (unsigned int k = 0; k < axes.z.bins; k++) {
                if (grid.index(i, j, k)) {
                    points[index] = Vector3<double>(axes.x.min + i*xstep, axes.y.min + j*ystep, axes.z.min + k*zstep);
                }
            }
        }
    }
    Fval::set_points(std::move(points));
}