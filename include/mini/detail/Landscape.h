#pragma once

#include <vector>

#include <dataset/SimpleDataset.h>
#include <plots/PlotOptions.h>
#include <mini/detail/Evaluation.h> // No forward declaration since the primary function of this class is to forward this object to the user

namespace mini {
    class Landscape : public plots::Plottable {
        public: 
            /**
             * @brief Default constructor.
             */
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