#pragma once

#include <vector>

#include <em/ImageStack.h>
#include <em/CullingStrategy.h>
#include <data/Protein.h>
#include <histogram/ScatteringHistogram.h>

namespace em {
    class PartialHistogramManager {
        public:
            PartialHistogramManager(const ImageStack& images);

            ~PartialHistogramManager() = default;

            histogram::ScatteringHistogram get_histogram(double cutoff);

            std::shared_ptr<Protein> get_protein() const;

            std::shared_ptr<Protein> get_protein(double cutoff);

            void set_cutoff_levels(std::vector<double> levels);

            /**
             * @brief Alternate slower approach to generating the histogram. 
             */
            histogram::ScatteringHistogram get_histogram_slow(double cutoff) const;
        
        private:
            const ImageStack& images; 
            std::shared_ptr<Protein> protein;
            std::vector<double> charge_levels = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 100000};
            double previous_cutoff = 0;

            void update_protein(double cutoff);

            std::unique_ptr<Protein>generate_protein(double cutoff) const;

            std::vector<Atom> generate_atoms(double cutoff) const;
    };
}