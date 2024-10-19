#pragma once

#include <data/DataFwd.h>

#include <vector>

namespace hist {
    /**
     * @brief Evaluate the exact Debye transform for the given molecule, ignoring waters, along the given q values.
     *
     * exp(-q*q) is used as the form factor.
     */
    std::vector<double> exact_debye_transform(const data::Molecule& molecule, const std::vector<double>& q_vals);
}