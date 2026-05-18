// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PairSchedule.h>

#include <array>
#include <cmath>
#include <map>

using namespace ausaxs;
using namespace ausaxs::symmetry;

namespace {
    // resolution at which two transforms are considered identical (1e-4)
    constexpr double key_resolution = 1e4;

    // canonical hashable key for a rigid transform: 9 rotation + 3 translation entries, rounded
    using TransformKey = std::array<long, 12>;

    TransformKey make_key(const Matrix<double>& R, const Vector3<double>& T) {
        TransformKey k;
        int idx = 0;
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
                k[idx++] = std::llround(R(i, j)*key_resolution);
            }
        }
        for (int i = 0; i < 3; ++i) {k[idx++] = std::llround(T[i]*key_resolution);}
        return k;
    }
}

std::vector<CopyPair> ausaxs::symmetry::compute_pair_schedule(const std::vector<AffineTransform>& placements) {
    int n = static_cast<int>(placements.size());
    std::map<TransformKey, CopyPair> buckets; // representative pair + running scale per class

    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            // relative transform D = placement_j^{-1} o placement_i; rigid, so inverse is the transpose
            Matrix<double> RjT = placements[j].rotation.transpose();
            Matrix<double> D_R = RjT*placements[i].rotation;
            Vector3<double> D_T = RjT*(placements[i].translation - placements[j].translation);

            // distance is symmetric: D and D^{-1} describe the same pair of bodies
            Matrix<double> Dinv_R = D_R.transpose();
            Vector3<double> Dinv_T = -(Dinv_R*D_T);

            TransformKey key = std::min(make_key(D_R, D_T), make_key(Dinv_R, Dinv_T));
            auto [it, inserted] = buckets.try_emplace(key, CopyPair{i, j, 1});
            if (!inserted) {++it->second.scale;}
        }
    }

    std::vector<CopyPair> schedule;
    schedule.reserve(buckets.size());
    for (const auto& [key, pair] : buckets) {schedule.push_back(pair);}
    return schedule;
}
