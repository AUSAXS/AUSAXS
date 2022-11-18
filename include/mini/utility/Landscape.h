#pragma once

#include <vector>
#include <list>
#include <string>

#include <mini/utility/Evaluation.h>
#include <math/Matrix.h>

namespace mini {
    /**
     * @brief A 2D grid of function evaluations. 
     *        The grid must be regularly spaced.
     */
    struct Landscape {
        /**
         * @brief Default constructor.
         */
        Landscape() noexcept {}

        /**
         * @brief Construct a Landscape from a saved file. 
         * 
         * @param file Path to the file. 
         */
        Landscape(std::string file);

        /**
         * @brief Construct an empty Landscape with a given number of dimensions.
         */
        Landscape(unsigned int xbins, unsigned int ybins) noexcept : x(xbins), y(ybins), z(xbins, ybins) {}

        /**
         * @brief Flip the x and y axes. 
         */
        void rotate90() noexcept;

        /**
         * @brief Add a set of evaluated points.
         */
        void add_evaluations(const std::vector<Evaluation>& evaluations) noexcept;

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

        std::list<Evaluation> evaluations;
        std::vector<double> x, y;
        Matrix<double> z;
    };    
}