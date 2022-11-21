#pragma once

#include <vector>

#include <mini/utility/Evaluation.h>
#include <utility/SimpleDataset.h>

namespace mini {
    class Landscape {
        public: 
            /**
             * @brief Default constructor.
             */
            Landscape() noexcept = default;

            Landscape(unsigned int size) : evals(size) {};

            Landscape(std::vector<Evaluation> evals) : evals(std::move(evals)) {}

            /**
             * @brief Convert this landscape to a SimpleDataset. 
             *        Requires that the landscape is 1-dimensional. 
             *        The dataset will be sorted by the x-axis. 
             */
            SimpleDataset as_dataset() const;

            void append(std::vector<Evaluation> evals) {this->evals.insert(this->evals.end(), evals.begin(), evals.end());}
            void append(Landscape evals) {append(evals.evals);}

            std::vector<Evaluation> evals;
    };    
}