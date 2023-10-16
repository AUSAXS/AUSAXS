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
#include <io/ExistingFile.h>
#include <constants/Constants.h>

#include <atomic>
#include <thread>
#include <random>
#include <fstream>
#include <csignal>
#include <mutex>

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

void CrystalScattering::random_rotation() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    double theta = std::acos(2*dis(gen) - 1);
    double phi = 2*M_PI*dis(gen);
    double psi = 2*M_PI*dis(gen);
    auto v = Vector3<double>(std::cos(theta)*std::cos(phi), std::sin(theta)*std::cos(phi), std::sin(phi));
    rotate(v, psi);
}

void CrystalScattering::rotate(const Vector3<double>& axis, double angle) {
    auto p = Fval::get_points();
    for (auto& point : p) {
        point.rotate(axis, angle);
    }
}

SimpleDataset CrystalScattering::rotational_average(unsigned int n) {
    std::vector<SimpleDataset> datasets;
    for (unsigned int i = 0; i < n; i++) {
        datasets.push_back(calculate());
        random_rotation();
    }

    SimpleDataset result;
    for (unsigned int i = 0; i < datasets[0].size_rows(); i++) {
        double sum = 0;
        for (auto& dataset : datasets) {
            sum += dataset.y(i);
        }
        result.push_back(datasets[0].x(i), sum/datasets.size());
    }
    return result;
}

CrystalScattering::CrystalScattering(const grid::Grid& grid) {
    initialize();
    convert_grid(grid);
}

void CrystalScattering::initialize() {
    miller_strategy = factory::construct_miller_strategy();
}

std::mutex checkpointMutex;
void CrystalScattering::save_checkpoint(unsigned int stop, const std::vector<Fval>& fvals) {
    std::lock_guard<std::mutex> lock(checkpointMutex); // Lock the mutex
    ::io::File file(settings::general::output + "temp/checkpoint.dat"); 
    file.create();

    std::ofstream checkpoint_file(file, std::ios::binary);
    if (!checkpoint_file.is_open()) {throw except::io_error("CrystalScattering::save_checkpoint: Could not open checkpoint file.");}

    // we need to know the current number of calculated points
    checkpoint_file.write(reinterpret_cast<const char*>(&stop), sizeof(stop));

    // we also save the total size of the fvals vector so we can make a simple consistency check when loading
    unsigned int fvals_size = fvals.size();
    checkpoint_file.write(reinterpret_cast<const char*>(&fvals_size), sizeof(fvals_size));

    // finally we write the first 'stop' elements of the fvals vector
    checkpoint_file.write(reinterpret_cast<const char*>(fvals.data()), stop*sizeof(Fval));
    checkpoint_file.close();
}

unsigned int CrystalScattering::load_checkpoint(std::vector<Fval>& fvals) {
    ::io::File file(settings::general::output + "temp/checkpoint.dat"); 

    std::ifstream checkpoint_file(file, std::ios::binary);
    if (!checkpoint_file.is_open()) {
        if (settings::general::verbose) {std::cout << "Could not load any points from checkpoint file.";}
        return 0;
    }

    // first we read the number of saved points
    unsigned int stop = 0;
    checkpoint_file.read(reinterpret_cast<char*>(&stop), sizeof(stop));
    if (stop > fvals.size()) {throw except::unexpected("CrystalScattering::load_checkpoint: incompatible checkpoint file. Did you change the settings?.");}

    // then we read the previous size of the fvals vector
    unsigned int fvals_size = 0;
    checkpoint_file.read(reinterpret_cast<char*>(&fvals_size), sizeof(fvals_size));
    if (fvals_size != fvals.size()) {throw except::unexpected("CrystalScattering::load_checkpoint: incompatible checkpoint file. Did you change the settings?.");}

    // finally we read the first 'stop' elements of the fvals vector
    checkpoint_file.read(reinterpret_cast<char*>(fvals.data()), stop*sizeof(Fval));
    checkpoint_file.close();

    if (settings::general::verbose) {std::cout << "Loaded " << stop << " points from checkpoint file." << std::endl;}
    return stop;
}

bool interrupt_signal = false;
void interrupt_handler(int signal) {
    std::cout << "Interrupt signal received. Finishing current calculations and writing a checkpoint before exiting. \nInterrupt again to exit immediately and lose current progress." << std::endl;
    std::signal(SIGINT, SIG_DFL);
    if (signal == SIGINT) {interrupt_signal = true;}
}

void interrupt_handler_save(int) {
    std::cout << "Program cannot be interrupted while writing a checkpoint file. Please wait." << std::endl;
}

SimpleDataset CrystalScattering::calculate() const {
    if (Fval::get_points().empty()) {throw except::invalid_argument("CrystalScattering::calculate: No points were set.");}
    if (Fval::get_basis().x.x() == 0 || Fval::get_basis().y.y() == 0 || Fval::get_basis().z.z() == 0) {throw except::invalid_argument("CrystalScattering::calculate: No basis was set.");}
    auto millers = miller_strategy->generate();

    std::vector<Fval> fvals(millers.size());
    std::atomic<unsigned int> index = load_checkpoint(fvals);

    auto dispatcher = [&] () {
        while (true) {
            unsigned int start = index.fetch_add(1000);
            unsigned int end = std::min<unsigned int>(start + 1000, millers.size());
            if (start >= millers.size() || interrupt_signal) {
                index = std::min(index.load(), end);
                break;
            }

            std::cout << start << "/" << millers.size() << "          \r" << std::flush;
            for (unsigned int i = start; i < end; i++) {
                fvals[i] = Fval(millers[i].h, millers[i].k, millers[i].l);
            }
        }
    };

    // if the checkpoint file does not contain all points, we need to calculate the remaining points
    if (index < millers.size()) {
        // register the interrupt signal handler. we need this to ensure that the checkpoint file is saved properly before the program exits
        std::signal(SIGINT, interrupt_handler);

        // start threads
        std::vector<std::thread> threads;
        for (unsigned int i = 0; i < settings::general::threads; i++) {
            threads.push_back(std::thread(dispatcher));
        }

        // wait for threads to finish
        for (auto& thread : threads) {
            thread.join();
        }

        // save final checkpoint
        std::signal(SIGINT, interrupt_handler_save);
        save_checkpoint(index, fvals);
        std::signal(SIGINT, SIG_DFL); // reset the interrupt signal handler
        if (interrupt_signal) {raise(SIGINT);} // raise the interrupt signal again to ensure that the program exits
    }

    // sort fvals by q
    std::sort(fvals.begin(), fvals.end(), [] (const Fval& a, const Fval& b) {return a.qlength < b.qlength;});

    // prepare the q profile
    std::vector<double> bins(constants::axes::q_axis.get_bin(settings::axes::qmax));
    double logmin = std::log10(settings::axes::qmin);
    double logmax = std::log10(settings::axes::qmax);
    double logstep = (logmax - logmin)/bins.size();
    for (unsigned int i = 0; i < bins.size(); i++) {
        bins[i] = std::pow(10, logmin + i*logstep);
    }

    // bin the data
    SimpleDataset data;
    unsigned int bin_index = 0;
    while (bin_index < fvals.size() && fvals[bin_index].qlength < bins[0]) {bin_index++;} // skip all values below the first bin
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