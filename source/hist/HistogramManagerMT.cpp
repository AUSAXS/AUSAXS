#include <data/Protein.h>
#include <hist/HistogramManagerMT.h>

#include <BS_thread_pool.hpp>

using namespace hist;

HistogramManagerMT::HistogramManagerMT(HistogramManager& hm) : HistogramManager(hm) {}

HistogramManagerMT::~HistogramManagerMT() = default;

Histogram HistogramManagerMT::calculate() {
    return HistogramManager::calculate().p;
}

ScatteringHistogram HistogramManagerMT::calculate_all() {
    auto atoms = protein->atoms();
    auto waters = protein->waters();

    double width = setting::axes::scattering_intensity_plot_binned_width;
    Axis axes = Axis(setting::axes::max_distance/width, 0, setting::axes::max_distance); 

    // create a more compact representation of the coordinates
    // extremely wasteful to calculate this from scratch every time
    std::vector<float> data_p(atoms.size()*4);
    for (unsigned int i = 0; i < atoms.size(); i++) {
        const Atom& a = atoms[i]; 
        data_p[4*i] = a.coords.x();
        data_p[4*i+1] = a.coords.y();
        data_p[4*i+2] = a.coords.z();
        data_p[4*i+3] = a.effective_charge*a.occupancy;
    }

    std::vector<float> data_h(waters.size()*4);
    for (unsigned int i = 0; i < waters.size(); i++) {
        const Water& a = waters[i]; 
        data_h[4*i] = a.coords.x();
        data_h[4*i+1] = a.coords.y();
        data_h[4*i+2] = a.coords.z();
        data_h[4*i+3] = a.effective_charge*a.occupancy;
    }

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    BS::thread_pool pool(setting::general::threads);
    auto calc_pp = [&data_p, &axes, &atoms, &width] (unsigned int imin, unsigned int imax) {
        std::vector<double> p_pp(axes.bins, 0);
        for (unsigned int i = imin; i < imax; i++) {
            for (unsigned int j = i+1; j < atoms.size(); j++) {
                float weight = data_p[4*i+3]*data_p[4*j+3]; // Z1*Z2*w1*w2
                float dx = data_p[4*i  ] - data_p[4*j  ];
                float dy = data_p[4*i+1] - data_p[4*j+1];
                float dz = data_p[4*i+2] - data_p[4*j+2];
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        return p_pp;
    };
    auto calc_hh = [&data_h, &axes, &waters, &width] (unsigned int imin, unsigned int imax) {
        std::vector<double> p_hh(axes.bins, 0);
        for (unsigned int i = imin; i < imax; i++) {
            for (unsigned int j = i+1; j < waters.size(); j++) {
                float weight = data_h[4*i+3]*data_h[4*j+3];
                float dx = data_h[4*i  ] - data_h[4*j  ];
                float dy = data_h[4*i+1] - data_h[4*j+1];
                float dz = data_h[4*i+2] - data_h[4*j+2];
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_hh[dist/width] += 2*weight;
            }
        }
        return p_hh;
    };
    auto calc_hp = [&data_h, &data_p, &axes, &waters, &atoms, &width] (unsigned int imin, unsigned int imax) {
        std::vector<double> p_hp(axes.bins, 0);
        for (unsigned int i = imin; i < imax; i++) {
            for (unsigned int j = 0; j < atoms.size(); j++) {
                float weight = data_h[4*i+3]*data_p[4*j+3];
                float dx = data_h[4*i  ] - data_p[4*j  ];
                float dy = data_h[4*i+1] - data_p[4*j+1];
                float dz = data_h[4*i+2] - data_p[4*j+2];
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_hp[dist/width] += 2*weight;
            }
        }
        return p_hp;
    };

    //##############//
    // SUBMIT TASKS //
    //##############//
    unsigned int job_size = setting::general::detail::job_size;
    BS::multi_future<std::vector<double>> pp;
    for (unsigned int i = 0; i < atoms.size(); i+=job_size) {
        pp.push_back(pool.submit(calc_pp, i, std::min(i+job_size, (unsigned int)atoms.size())));
    }
    BS::multi_future<std::vector<double>> hh;
    for (unsigned int i = 0; i < waters.size(); i+=job_size) {
        hh.push_back(pool.submit(calc_hh, i, std::min(i+job_size, (unsigned int)waters.size())));
    }
    BS::multi_future<std::vector<double>> hp;
    for (unsigned int i = 0; i < waters.size(); i+=job_size) {
        hp.push_back(pool.submit(calc_hp, i, std::min(i+job_size, (unsigned int)waters.size())));
    }
    utility::print_info("Submitted " + std::to_string(pp.size() + hh.size() + hp.size()) + " tasks to thread pool. Waiting for completion...");
    pool.wait_for_tasks();

    //#################//
    // COLLECT RESULTS //
    //#################//    
    std::vector<double> p_pp(axes.bins, 0);
    for (const auto& tmp : pp.get()) {
        for (unsigned int j = 0; j < axes.bins; j++) {
            p_pp[j] += tmp[j];
        }
    }
    std::vector<double> p_hh(axes.bins, 0);
    for (const auto& tmp : hh.get()) {
        for (unsigned int j = 0; j < axes.bins; j++) {
            p_hh[j] += tmp[j];
        }
    }
    std::vector<double> p_hp(axes.bins, 0);
    for (const auto& tmp : hp.get()) {
        for (unsigned int j = 0; j < axes.bins; j++) {
            p_hp[j] += tmp[j];
        }
    }

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    for (unsigned int i = 0; i < atoms.size(); i++) {p_pp[0] += data_p[4*i+3]*data_p[4*i+3];}
    for (unsigned int i = 0; i < waters.size(); i++) {p_hh[0] += data_h[4*i+3]*data_h[4*i+3];}

    // downsize our axes to only the relevant area
    int max_bin = 10; // minimum size is 10
    for (int i = axes.bins-1; i >= 10; i--) {
        if (p_pp[i] != 0 || p_hh[i] != 0 || p_hp[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }
    p_pp.resize(max_bin);
    p_hh.resize(max_bin);
    p_hp.resize(max_bin);
    axes = Axis{max_bin, 0, max_bin*width}; 

    // calculate p_tot    
    std::vector<double> p_tot(max_bin, 0);
    for (int i = 0; i < max_bin; i++) {p_tot[i] = p_pp[i] + p_hh[i] + p_hp[i];}

    return ScatteringHistogram(p_pp, p_hh, p_hp, p_tot, axes);
}
