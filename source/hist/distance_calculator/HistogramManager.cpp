#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/state/StateManager.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/detail/CompactCoordinates.h>
#include <settings/HistogramSettings.h>
#include <constants/Constants.h>

using namespace hist;

HistogramManager::HistogramManager(view_ptr<const data::Molecule> protein) : BodyTracker(protein), protein(protein) {}

HistogramManager::HistogramManager(const HistogramManager& hm) : BodyTracker(hm.protein), protein(hm.protein) {}

HistogramManager::~HistogramManager() = default;

std::unique_ptr<DistanceHistogram> HistogramManager::calculate() {return calculate_all();}

std::unique_ptr<CompositeDistanceHistogram> HistogramManager::calculate_all() {
    std::vector<double> p_pp( constants::axes::d_axis.bins, 0);
    std::vector<double> p_hh( constants::axes::d_axis.bins, 0);
    std::vector<double> p_hp( constants::axes::d_axis.bins, 0);

    hist::detail::CompactCoordinates data_p(protein->get_bodies());
    hist::detail::CompactCoordinates data_h = hist::detail::CompactCoordinates(protein->get_waters());

    // calculate p-p distances
    for (unsigned int i = 0; i < data_p.get_size(); ++i) {
        unsigned int j = i+1;
        for (; j+7 < data_p.get_size(); j+=8) {
            auto res = data_p[i].evaluate(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3], data_p[j+4], data_p[j+5], data_p[j+6], data_p[j+7]);
            for (unsigned int k = 0; k < 8; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j+3 < data_p.get_size(); j+=4) {
            auto res = data_p[i].evaluate(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3]);
            for (unsigned int k = 0; k < 4; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j < data_p.get_size(); ++j) {
            auto res = data_p[i].evaluate(data_p[j]);
            p_pp[res.distance] += 2*res.weight;
        }
    }

    // add self-correlation
    p_pp[0] = std::accumulate(data_p.get_data().begin(), data_p.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + std::pow(val.value.w, 2);} );

    for (unsigned int i = 0; i < data_h.get_size(); ++i) {
        // calculate h-h distances
        {
            unsigned int j = i+1;
            for (; j+7 < data_h.get_size(); j+=8) {
                auto res = data_h[i].evaluate(data_h[j], data_h[j+1], data_h[j+2], data_h[j+3], data_h[j+4], data_h[j+5], data_h[j+6], data_h[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {p_hh[res.distance[k]] += 2*res.weight[k];}
            }

            for (; j+3 < data_h.get_size(); j+=4) {
                auto res = data_h[i].evaluate(data_h[j], data_h[j+1], data_h[j+2], data_h[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {p_hh[res.distance[k]] += 2*res.weight[k];}
            }

            for (; j < data_h.get_size(); ++j) {
                auto res = data_h[i].evaluate(data_h[j]);
                p_hh[res.distance] += 2*res.weight;
            }
        }
        
        // calculate h-p distances
        {
            unsigned int j = 0;
            for (; j+7 < data_p.get_size(); j+=8) {
                auto res = data_h[i].evaluate(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3], data_p[j+4], data_p[j+5], data_p[j+6], data_p[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {p_hp[res.distance[k]] += res.weight[k];}
            }

            for (; j+3 < data_p.get_size(); j+=4) {
                auto res = data_h[i].evaluate(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {p_hp[res.distance[k]] += res.weight[k];}
            }

            for (; j < data_p.get_size(); ++j) {
                auto res = data_h[i].evaluate(data_p[j]);
                p_hp[res.distance] += res.weight;
            }
        }
    }

    // add self-correlation
    p_hh[0] = std::accumulate(data_h.get_data().begin(), data_h.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + std::pow(val.value.w, 2);} );

    // calculate p_tot
    std::vector<double> p_tot(constants::axes::d_axis.bins, 0);
    for (unsigned int i = 0; i < p_pp.size(); ++i) {p_tot[i] = p_pp[i] + p_hh[i] + 2*p_hp[i];}

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (int i = p_tot.size()-1; i >= 10; i--) {
        if (p_tot[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p_pp.resize(max_bin);
    p_hh.resize(max_bin);
    p_hp.resize(max_bin);
    p_tot.resize(max_bin);
    return std::make_unique<CompositeDistanceHistogram>(
        std::move(p_pp), 
        std::move(p_hp), 
        std::move(p_hh), 
        std::move(p_tot), 
        Axis(0, max_bin*constants::axes::d_axis.width(), max_bin)
    );
}