#pragma once

#include <mini/detail/Landscape.h>
#include <mini/detail/Evaluation.h>
#include <math/Matrix.h>

#include <string>

namespace io {class ExistingFile;}
namespace mini {
    struct RegularLandscape : Landscape {
        RegularLandscape(const Landscape& l);

        /**
         * @brief Construct a Landscape from a saved file. 
         * 
         * @param file Path to the file. 
         */
        RegularLandscape(const io::ExistingFile& file);

        /**
         * @brief Flip the x and y axes. 
         */
        void rotate90() noexcept;

        /**
         * @brief Find the minimum evaluated point in the landscape.
         */
        [[nodiscard]] Evaluation find_min_eval() const;

        /**
         * @brief Find the minimum value in the landscape.
         */
        [[nodiscard]] Evaluation find_min_val() const;

        /**
         * @brief Save the landscape to a file.
         *        Only the raw data is saved, not the evaluated points.
         */
        void save(std::string filename) const;

        /**
         * @brief Load a landscape from a file.
         * 
         * @param filename Path to the file.
         */
        void load(std::string filename);

        std::vector<double> x, y;
        Matrix<double> z;
    };    
}