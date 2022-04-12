#pragma once

#include <vector>

#include <em/ImageStack.h>
#include <em/CullingStrategy.h>
#include <data/Protein.h>
#include <ScatteringHistogram.h>

namespace em {
    class PartialHistogramManager {
        public:
            PartialHistogramManager(const ImageStack& images);

            ScatteringHistogram get_histogram(double cutoff);

            std::shared_ptr<Protein> get_protein() const;

            void set_cutoff_levels(std::vector<double> levels);

            void update_protein(double cutoff);

            /**
             * @brief Alternate slower approach to generating the histogram. 
             */
            ScatteringHistogram get_histogram_slow(double cutoff) const;
        
        private:
            const ImageStack& images; 
            std::unique_ptr<CullingStrategy> culler;
            std::shared_ptr<Protein> protein;
            std::vector<double> charge_levels = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 100000};
            double previous_cutoff = 0;

            std::unique_ptr<Protein>generate_protein(double cutoff) const;

            std::vector<Atom> generate_atoms(double cutoff) const;
    };
}