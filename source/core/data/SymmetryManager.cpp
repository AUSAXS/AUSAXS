// #include <data/SymmetryManager.h>
// #include <data/Molecule.h>
// #include <data/Body.h>
// #include <data/record/Atom.h>
// #include <data/record/Water.h>
// #include <hist/distance_calculator/SimpleCalculator.h>
// #include <hist/distribution/GenericDistribution1D.h>
// #include <hist/intensity_calculator/CompositeDistanceHistogram.h>

// #include <cassert>

// using namespace ausaxs;

// enum type {
//     SELF_AA, SELF_AW, SELF_WW, 
//     SELF_SYM_AA, SELF_SYM_AW, SELF_SYM_WW,
//     CROSS_AA, CROSS_AW, CROSS_WW
// };

// template<bool use_weighted_distribution>
// std::unique_ptr<hist::CompositeDistanceHistogram> hist::detail::SymmetryManager::calculate(
//     const data::Molecule& protein
// ) {
//     using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
//     hist::distance_calculator::SimpleCalculator<use_weighted_distribution> calculator;

//     std::vector<type> self_indices, cross_indices;
//     std::vector<int> needs_scaling;
//     for (unsigned int i = 0; i < protein.size_body(); ++i) {
//         const auto& body = protein.get_body(i);
//         const auto& atoms = body.get_atoms();
//         const auto& waters = body.get_waters();

//         hist::detail::CompactCoordinates data_a(atoms);
//         hist::detail::CompactCoordinates data_w(waters);

//         // base histogram
//         calculator.enqueue_calculate_self(data_a);
//         calculator.enqueue_calculate_self(data_w);
//         needs_scaling.push_back(calculator.enqueue_calculate_cross(data_a, data_w));
//         self_indices.push_back(SELF_AA);
//         self_indices.push_back(SELF_WW);
//         cross_indices.push_back(SELF_AW);

//         Vector3<double> cm = body.get_cm();
//         for (const auto& symmetry : body.get_symmetries()) {

//             // internal histograms
//             Matrix<double> rotate = matrix::rotation_matrix(symmetry.rotate.x(), symmetry.rotate.y(), symmetry.rotate.z());
//             auto transform = [&symmetry, &rotate, &cm] (const auto& a) {
//                 return rotate*(a.get_coordinates()-cm) + cm + symmetry.translate;
//             };

//             std::vector<data::record::Atom> atoms_transformed(atoms.size());
//             std::vector<data::record::Water> waters_transformed(waters.size());
//             std::transform(atoms.begin(), atoms.end(), atoms_transformed.begin(), transform);
//             std::transform(waters.begin(), waters.end(), waters_transformed.begin(), transform);

//             hist::detail::CompactCoordinates data_a_transformed(atoms_transformed);
//             hist::detail::CompactCoordinates data_w_transformed(waters_transformed);

//             calculator.enqueue_calculate_cross(data_a, data_a_transformed);
//             calculator.enqueue_calculate_cross(data_a, data_w_transformed);
//             cross_indices.push_back(SELF_SYM_AA);
//             cross_indices.push_back(SELF_SYM_AW);

//             calculator.enqueue_calculate_cross(data_w, data_w_transformed);
//             calculator.enqueue_calculate_cross(data_w, data_a_transformed);
//             cross_indices.push_back(SELF_SYM_WW);
//             cross_indices.push_back(SELF_SYM_AW);

//             // external histograms with other bodies
//             for (unsigned int j = 0; j < protein.size_body(); ++j) {
//                 if (j == i) {continue;}

//                 const auto& body2 = protein.get_body(j);
//                 const auto& atoms2 = body2.get_atoms();
//                 const auto& waters2 = body2.get_waters();

//                 hist::detail::CompactCoordinates data_a2(atoms2);
//                 hist::detail::CompactCoordinates data_w2(waters2);

//                 calculator.enqueue_calculate_cross(data_a2, data_a_transformed);
//                 calculator.enqueue_calculate_cross(data_a2, data_w_transformed);
//                 cross_indices.push_back(CROSS_AA);
//                 cross_indices.push_back(CROSS_AW);

//                 calculator.enqueue_calculate_cross(data_w2, data_w_transformed);
//                 calculator.enqueue_calculate_cross(data_w2, data_a_transformed);
//                 cross_indices.push_back(CROSS_WW);
//                 cross_indices.push_back(CROSS_AW);
//             }
//         }

//         // external histograms with other bodies
//         for (unsigned int j = 0; j < protein.size_body(); ++j) {
//             if (j == i) {continue;}

//             const auto& body2 = protein.get_body(j);
//             auto atoms2 = body2.get_atoms();
//             auto waters2 = body2.get_waters();

//             hist::detail::CompactCoordinates data_a2(atoms2);
//             hist::detail::CompactCoordinates data_w2(waters2);

//             calculator.enqueue_calculate_cross(data_a, data_a2);
//             calculator.enqueue_calculate_cross(data_a, data_w2);
//             cross_indices.push_back(CROSS_AA);
//             cross_indices.push_back(CROSS_AW);

//             calculator.enqueue_calculate_cross(data_w, data_w2);
//             calculator.enqueue_calculate_cross(data_w, data_a2);
//             cross_indices.push_back(CROSS_WW);
//             cross_indices.push_back(CROSS_AW);
//         }
//     }

//     auto res = calculator.run();

//     assert(res.self.size() == self_indices.size() && "SymmetryManager::calculate: self size mismatch");
//     assert(res.cross.size() == cross_indices.size() && "SymmetryManager::calculate: cross size mismatch");

//     GenericDistribution1D_t p_aa = res.self[0]*protein.get_body(0).size_symmetry();
//     GenericDistribution1D_t p_ww = res.self[1]*protein.get_body(0).size_symmetry();
//     GenericDistribution1D_t p_aw = res.cross[0]*protein.get_body(0).size_symmetry();
//     for (int i = 1; i < static_cast<int>(protein.size_body()); ++i) {
//         int duplicates = protein.get_body(i).size_symmetry();
//         p_aa += res.self[2*i]*duplicates;
//         p_ww += res.self[2*i+1]*duplicates;
//         res.cross[needs_scaling[i]] *= duplicates;
//     }

//     for (int i = 1; i < static_cast<int>(cross_indices.size()); ++i) {
//         switch (cross_indices[i]) {
//             case SELF_AA:
//             case CROSS_AA:
//             case SELF_SYM_AA: p_aa += res.cross[i]; break;

//             case SELF_WW:
//             case CROSS_WW:
//             case SELF_SYM_WW: p_ww += res.cross[i]; break;

//             case SELF_AW:
//             case CROSS_AW:
//             case SELF_SYM_AW: p_aw += res.cross[i]; break;
//         }
//     }

//     // calculate p_tot
//     GenericDistribution1D_t p_tot(constants::axes::d_axis.bins);
//     for (unsigned int i = 0; i < p_tot.size(); ++i) {p_tot.index(i) = p_aa.index(i) + p_ww.index(i) + 2*p_aw.index(i);}

//     // downsize our axes to only the relevant area
//     unsigned int max_bin = 10; // minimum size is 10
//     for (int i = p_tot.size()-1; i >= 10; i--) {
//         if (p_tot.index(i) != 0) {
//             max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
//             break;
//         }
//     }
//     p_aa.resize(max_bin);
//     p_ww.resize(max_bin);
//     p_aw.resize(max_bin);
//     p_tot.resize(max_bin);

//     if constexpr (use_weighted_distribution) {
//         return std::make_unique<CompositeDistanceHistogram>(
//             Distribution1D(std::move(p_aa)), 
//             Distribution1D(std::move(p_aw)), 
//             Distribution1D(std::move(p_ww)), 
//             std::move(p_tot)
//         );
//     } else {
//         return std::make_unique<CompositeDistanceHistogram>(
//             std::move(p_aa), 
//             std::move(p_aw), 
//             std::move(p_ww), 
//             std::move(p_tot)
//         );
//     }
// }

// template std::unique_ptr<hist::CompositeDistanceHistogram> hist::detail::SymmetryManager::calculate<true>(const data::Molecule&);
// template std::unique_ptr<hist::CompositeDistanceHistogram> hist::detail::SymmetryManager::calculate<false>(const data::Molecule&);