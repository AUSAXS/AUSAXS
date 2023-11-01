#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/state/StateManager.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <settings/HistogramSettings.h>
#include <constants/Constants.h>

using namespace hist;

template<bool use_weighted_distribution>
HistogramManager<use_weighted_distribution>::HistogramManager(view_ptr<const data::Molecule> protein) : IHistogramManager(protein), protein(protein) {}

template<bool use_weighted_distribution>
HistogramManager<use_weighted_distribution>::HistogramManager(const HistogramManager& hm) : IHistogramManager(hm.protein), protein(hm.protein) {}

template<bool use_weighted_distribution>
HistogramManager<use_weighted_distribution>::~HistogramManager() = default;

template<bool use_weighted_distribution>
std::unique_ptr<DistanceHistogram> HistogramManager<use_weighted_distribution>::calculate() {return calculate_all();}

inline void add8_impl(hist::Distribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d, const hist::detail::CompactCoordinatesData& e, const hist::detail::CompactCoordinatesData& f, const hist::detail::CompactCoordinatesData& g, const hist::detail::CompactCoordinatesData& h) {
    auto res = data.evaluate(a, b, c, d, e, f, g, h);
    for (unsigned int k = 0; k < 8; ++k) {p.add(res.distance[k], res.weight[k]);}
}

inline void add8_impl(hist::WeightedDistribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d, const hist::detail::CompactCoordinatesData& e, const hist::detail::CompactCoordinatesData& f, const hist::detail::CompactCoordinatesData& g, const hist::detail::CompactCoordinatesData& h) {
    auto res = data.evaluate_rounded(a, b, c, d, e, f, g, h);
    for (unsigned int k = 0; k < 8; ++k) {p.add(res.distance[k], res.weight[k]);}
}

template<bool use_weighted_distribution>
inline void add8(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d, const hist::detail::CompactCoordinatesData& e, const hist::detail::CompactCoordinatesData& f, const hist::detail::CompactCoordinatesData& g, const hist::detail::CompactCoordinatesData& h) {
    add8_impl(p, data, a, b, c, d, e, f, g, h);
}

inline void add4_impl(hist::Distribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d) {
    auto res = data.evaluate(a, b, c, d);
    for (unsigned int k = 0; k < 4; ++k) {p.add(res.distance[k], res.weight[k]);}
}

inline void add4_impl(hist::WeightedDistribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d) {
    auto res = data.evaluate_rounded(a, b, c, d);
    for (unsigned int k = 0; k < 4; ++k) {p.add(res.distance[k], res.weight[k]);}
}

template<bool use_weighted_distribution>
inline void add4(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d) {
    add4_impl(p, data, a, b, c, d);
}

inline void add1_impl(hist::WeightedDistribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a) {
    auto res = data.evaluate(a);
    p.add(res.distance, res.weight);
}

inline void add1_impl(hist::Distribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a) {
    auto res = data.evaluate_rounded(a);
    p.add(res.distance, res.weight);
}

template<bool use_weighted_distribution>
inline void add1(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a) {
    add1_impl(p, data, a);
}

template<bool use_weighted_distribution>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManager<use_weighted_distribution>::calculate_all() {
    typename hist::GenericDistribution1D<use_weighted_distribution>::type p_pp(constants::axes::d_axis.bins, 0);
    typename hist::GenericDistribution1D<use_weighted_distribution>::type p_hh(constants::axes::d_axis.bins, 0);
    typename hist::GenericDistribution1D<use_weighted_distribution>::type p_hp(constants::axes::d_axis.bins, 0);

    hist::detail::CompactCoordinates data_p(protein->get_bodies());
    hist::detail::CompactCoordinates data_h = hist::detail::CompactCoordinates(protein->get_waters());

    // calculate p-p distances
    for (unsigned int i = 0; i < data_p.get_size(); ++i) {
        unsigned int j = i+1;
        for (; j+7 < data_p.get_size(); j+=8) {
            add8<use_weighted_distribution>(p_pp, data_p[i], data_p[j], data_p[j+1], data_p[j+2], data_p[j+3], data_p[j+4], data_p[j+5], data_p[j+6], data_p[j+7]);
            // auto res = data_p[i].evaluate_rounded(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3], data_p[j+4], data_p[j+5], data_p[j+6], data_p[j+7]);
            // for (unsigned int k = 0; k < 8; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j+3 < data_p.get_size(); j+=4) {
            add4<use_weighted_distribution>(p_pp, data_p[i], data_p[j], data_p[j+1], data_p[j+2], data_p[j+3]);
            // auto res = data_p[i].evaluate_rounded(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3]);
            // for (unsigned int k = 0; k < 4; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j < data_p.get_size(); ++j) {
            add1<use_weighted_distribution>(p_pp, data_p[i], data_p[j]);
            // auto res = data_p[i].evaluate_rounded(data_p[j]);
            // p_pp[res.distance] += 2*res.weight;
        }
    }

    for (unsigned int i = 0; i < data_h.get_size(); ++i) {
        // calculate h-h distances
        {
            unsigned int j = i+1;
            for (; j+7 < data_h.get_size(); j+=8) {
                add8<use_weighted_distribution>(p_hh, data_h[i], data_h[j], data_h[j+1], data_h[j+2], data_h[j+3], data_h[j+4], data_h[j+5], data_h[j+6], data_h[j+7]);
                // auto res = data_h[i].evaluate_rounded(data_h[j], data_h[j+1], data_h[j+2], data_h[j+3], data_h[j+4], data_h[j+5], data_h[j+6], data_h[j+7]);
                // for (unsigned int k = 0; k < 8; ++k) {p_hh[res.distance[k]] += 2*res.weight[k];}
            }

            for (; j+3 < data_h.get_size(); j+=4) {
                add4<use_weighted_distribution>(p_hh, data_h[i], data_h[j], data_h[j+1], data_h[j+2], data_h[j+3]);
                // auto res = data_h[i].evaluate_rounded(data_h[j], data_h[j+1], data_h[j+2], data_h[j+3]);
                // for (unsigned int k = 0; k < 4; ++k) {p_hh[res.distance[k]] += 2*res.weight[k];}
            }

            for (; j < data_h.get_size(); ++j) {
                add1<use_weighted_distribution>(p_hh, data_h[i], data_h[j]);
                // auto res = data_h[i].evaluate_rounded(data_h[j]);
                // p_hh[res.distance] += 2*res.weight;
            }
        }
        
        // calculate h-p distances
        {
            unsigned int j = 0;
            for (; j+7 < data_p.get_size(); j+=8) {
                add8<use_weighted_distribution>(p_hp, data_h[i], data_p[j], data_p[j+1], data_p[j+2], data_p[j+3], data_p[j+4], data_p[j+5], data_p[j+6], data_p[j+7]);
                // auto res = data_h[i].evaluate_rounded(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3], data_p[j+4], data_p[j+5], data_p[j+6], data_p[j+7]);
                // for (unsigned int k = 0; k < 8; ++k) {p_hp[res.distance[k]] += res.weight[k];}
            }

            for (; j+3 < data_p.get_size(); j+=4) {
                add4<use_weighted_distribution>(p_hp, data_h[i], data_p[j], data_p[j+1], data_p[j+2], data_p[j+3]);
                // auto res = data_h[i].evaluate_rounded(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3]);
                // for (unsigned int k = 0; k < 4; ++k) {p_hp[res.distance[k]] += res.weight[k];}
            }

            for (; j < data_p.get_size(); ++j) {
                add1<use_weighted_distribution>(p_hp, data_h[i], data_p[j]);
                // auto res = data_h[i].evaluate_rounded(data_p[j]);
                // p_hp[res.distance] += res.weight;
            }
        }
    }

    // add self-correlation
    p_pp.add(0, std::accumulate(data_p.get_data().begin(), data_p.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + std::pow(val.value.w, 2);}));
    p_hh.add(0, std::accumulate(data_h.get_data().begin(), data_h.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + std::pow(val.value.w, 2);}));

    // calculate p_tot
    Distribution1D p_tot(constants::axes::d_axis.bins, 0);
    for (unsigned int i = 0; i < p_pp.size(); ++i) {p_tot.index(i) = p_pp.index(i) + p_hh.index(i) + 2*p_hp.index(i);}

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (int i = p_tot.size()-1; i >= 10; i--) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p_pp.resize(max_bin);
    p_hh.resize(max_bin);
    p_hp.resize(max_bin);
    p_tot.resize(max_bin);
    return std::make_unique<CompositeDistanceHistogram<use_weighted_distribution>>(
        std::move(p_pp), 
        std::move(p_hp), 
        std::move(p_hh), 
        std::move(p_tot), 
        Axis(0, max_bin*constants::axes::d_axis.width(), max_bin)
    );
}

template class hist::HistogramManager<false>;
template class hist::HistogramManager<true>;