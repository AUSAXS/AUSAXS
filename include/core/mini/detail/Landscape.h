// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <dataset/SimpleDataset.h>
#include <mini/detail/Evaluation.h>

#include <vector>

namespace ausaxs::mini {
    class Landscape {
        public: 
            Landscape() noexcept = default;

            Landscape(unsigned int size) : evals(size) {}

            Landscape(const std::vector<Evaluation>& evals);

            /**
             * @brief Convert this landscape to a SimpleDataset. 
             *        Requires that the landscape is 1-dimensional. 
             *        The dataset will be sorted by the x-axis. 
             */
            SimpleDataset as_dataset() const;

            void append(const std::vector<Evaluation>& evals);
            void append(const Landscape& evals);

            std::string to_string() const;

            std::vector<Evaluation> evals;
    };    
}