#pragma once

#include <vector>

#include <em/ImageStack.h>
#include <em/CullingStrategy.h>
#include <data/Protein.h>
#include <histogram/ScatteringHistogram.h>

namespace em {
    /**
     * @brief A helper class for the ImageStack. 
     * 
     * This class generates and updates histograms in a smart way. This is done by splitting the 
     * generated atoms into multiple bodies, and then utilizing the smart histogram manager from the Protein. 
     */
    class PartialHistogramManager {
        public:
            /**
             * @brief Construct a Manager from an ImageStack.
             */
            PartialHistogramManager(const ImageStack& images);

            /**
             * @brief Destructor.
             */
            ~PartialHistogramManager() = default;

            /**
             * @brief Get the histogram for a given cutoff.
             */
            hist::ScatteringHistogram get_histogram(double cutoff);

            /**
             * @brief Get the Protein backing this object. 
             */
            std::shared_ptr<Protein> get_protein() const;

            /**
             * @brief Get the Protein generated from a given cutoff.
             */
            std::shared_ptr<Protein> get_protein(double cutoff);

            /**
             * @brief Set the cutoff levels. 
             */
            void set_cutoff_levels(std::vector<double> levels);

            /**
             * @brief Alternate slower approach to generating the histogram. 
             *        This is primarily meant for debugging purposes. 
             */
            hist::ScatteringHistogram get_histogram_slow(double cutoff) const;
        
        private:
            const ImageStack& images; 
            std::shared_ptr<Protein> protein;
            std::vector<double> charge_levels = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 100000};
            double previous_cutoff = 0;

            /**
             * @brief Update the Protein to reflect a new cutoff value.
             */
            void update_protein(double cutoff);

            /**
             * @brief Generate a new Protein for a given cutoff. 
             */
            std::unique_ptr<Protein>generate_protein(double cutoff) const;

            /**
             * @brief Generate the atmos for a given cutoff.
             */
            std::vector<Atom> generate_atoms(double cutoff) const;
    };
}