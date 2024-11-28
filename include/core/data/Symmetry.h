#pragma once

#include <math/Matrix.h>

#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>

namespace ausaxs::hist::detail {
    class Symmetry {
        public:
            template<bool use_weighted_distribution>
            void calculate(
                hist::distance_calculator::SimpleCalculator<use_weighted_distribution>& calculator,
                const data::Molecule& protein
            ) {
                using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
                auto atoms = protein.get_atoms();
                auto waters = protein.get_waters();

                hist::detail::CompactCoordinates data_a(atoms);
                hist::detail::CompactCoordinates data_w(protein.get_waters());

                calculator.enqueue_calculate_self(data_a);
                calculator.enqueue_calculate_self(data_w);
                calculator.enqueue_calculate_cross(data_a, data_w);

                for (auto& transform : symmetries) {
                    std::vector<data::record::Atom> atoms_transformed(atoms.size());
                    std::vector<data::record::Water> waters_transformed(waters.size());
                    std::transform(atoms.begin(), atoms.end(), atoms_transformed.begin(), [&transform] (const data::record::Atom& a) {return transform*a.coords;});
                    std::transform(waters.begin(), waters.end(), waters_transformed.begin(), [&transform] (const data::record::Water& w) {return transform*w.coords;});

                    hist::detail::CompactCoordinates data_a_transformed(atoms_transformed);
                    hist::detail::CompactCoordinates data_w_transformed(waters_transformed);

                    calculator.enqueue_calculate_cross(data_a, data_a_transformed);
                    calculator.enqueue_calculate_cross(data_w, data_w_transformed);

                    calculator.enqueue_calculate_cross(data_a, data_w_transformed);
                    calculator.enqueue_calculate_cross(data_w, data_a_transformed);
                }
                auto res = calculator.run();

                // scale self-correlation
                GenericDistribution1D_t p_aa = res.self[0]*symmetries.size();
                GenericDistribution1D_t p_ww = res.self[1]*symmetries.size();
                GenericDistribution1D_t p_aw = res.cross[0];

                // scale cross-correlation
                for (unsigned int i = 1; i < res.cross.size(); ++i) {
                    p_aw += res.cross[i];
                }

                // calculate p_tot
                GenericDistribution1D_t p_tot(constants::axes::d_axis.bins);
                for (unsigned int i = 0; i < p_tot.size(); ++i) {p_tot.index(i) = p_aa.index(i) + p_ww.index(i) + 2*p_aw.index(i);}

                // downsize our axes to only the relevant area
                unsigned int max_bin = 10; // minimum size is 10
                for (int i = p_tot.size()-1; i >= 10; i--) {
                    if (p_tot.index(i) != 0) {
                        max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
                        break;
                    }
                }
                p_aa.resize(max_bin);
                p_ww.resize(max_bin);
                p_aw.resize(max_bin);
                p_tot.resize(max_bin);

                if constexpr (use_weighted_distribution) {
                    return std::make_unique<CompositeDistanceHistogram>(
                        Distribution1D(std::move(p_aa)), 
                        Distribution1D(std::move(p_aw)), 
                        Distribution1D(std::move(p_ww)), 
                        std::move(p_tot)
                    );
                } else {
                    return std::make_unique<CompositeDistanceHistogram>(
                        std::move(p_aa), 
                        std::move(p_aw), 
                        std::move(p_ww), 
                        std::move(p_tot)
                    );
                }
            }

        private:
            std::vector<Matrix<double>> symmetries;
    };
}