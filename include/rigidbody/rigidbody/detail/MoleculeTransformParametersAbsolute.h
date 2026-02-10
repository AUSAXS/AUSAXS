// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::detail {
    /**
     * @brief A small structure for storing the absolute transform parameters of a molecule. 
     */
    struct MoleculeTransformParametersAbsolute {
        MoleculeTransformParametersAbsolute();
        MoleculeTransformParametersAbsolute(observer_ptr<const Rigidbody> rigidbody) noexcept;
        MoleculeTransformParametersAbsolute(observer_ptr<const Rigidbody> rigidbody, double chi2) noexcept;
        ~MoleculeTransformParametersAbsolute();

        std::vector<parameter::BodyTransformParametersAbsolute> parameters;
        std::vector<data::Water> waters;
        double chi2 = 1e9;
    };
}