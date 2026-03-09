// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/ExactDebyeCalculator.h>
#include <data/Molecule.h>
#include <hist/detail/CompactCoordinates.h>

#include <cmath>
#include <span>

using namespace ausaxs;

std::vector<double> hist::exact_debye_transform(const data::Molecule& molecule, const std::vector<double>& q_vals) {
    using CC = hist::detail::CompactCoordinates<false>;
    auto data = CC(molecule.get_bodies());
    using ElemType = std::remove_reference_t<decltype(data[0])>;

    auto contribution = [] (double qr, float w) -> double {
        if (qr < 1e-9) {
            return w;
        } else {
            return w*std::sin(qr)/qr;
        }
    };

    std::vector<double> I;
    I.reserve(q_vals.size());
    for (const auto& q : q_vals) {
        double sum = 0;
        for (unsigned int i = 0; i < data.size(); ++i) {
            unsigned int j = i+1;
            for (; j+7 < data.size(); j+=8) {
                auto res = data[i].evaluate_8(std::span<const ElemType, 8>(&data[j], 8));
                for (unsigned int k = 0; k < 8; ++k) {
                    sum += contribution(q*res.distances[k], 2*res.weights[k]);
                }
            }

            for (; j+3 < data.size(); j+=4) {
                auto res = data[i].evaluate_4(std::span<const ElemType, 4>(&data[j], 4));
                for (unsigned int k = 0; k < 4; ++k) {
                    sum += contribution(q*res.distances[k], 2*res.weights[k]);
                }
            }

            for (; j < data.size(); ++j) {
                auto res = data[i].evaluate(data[j]);
                    sum += contribution(q*res.distance, 2*res.weight);
            }

            sum += std::pow(data[i].value.w, 2);
        }

        sum *= std::exp(-q*q);
        I.push_back(sum);
    }
    return I;
}