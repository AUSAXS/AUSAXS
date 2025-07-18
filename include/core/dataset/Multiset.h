// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <dataset/Dataset2D.h>
#include <io/ExistingFile.h>

#include <unordered_map>
#include <vector>
#include <string>

namespace ausaxs {
    class Multiset {
        public:
            [[deprecated]] Multiset() = default;

            [[deprecated]] Multiset(const io::Folder& path) {read(path);}

            [[deprecated]] explicit Multiset(unsigned int size) : data(size) {}

            [[deprecated]] explicit Multiset(const std::vector<Dataset2D>& data);

            [[deprecated]] explicit Multiset(const std::vector<SimpleDataset>& data);

            [[deprecated]] explicit Multiset(const Dataset2D& data);

            [[deprecated]] Multiset(const Dataset2D& data1, const Dataset2D& data2);

            const Dataset2D& operator[](unsigned int i) const;
            Dataset2D& operator[](unsigned int i);

            const Dataset2D& get_data(const std::string& name) const;
            Dataset2D& get_data(const std::string& name);

            const Dataset2D& get_data(unsigned int i) const;
            Dataset2D& get_data(unsigned int i);

            /**
             * @brief Get the number of Datasets contained in this Multiset. 
             */
            std::size_t size() const;

            /**
             * @brief Check if this Multifram is empty.
             */
            bool empty() const;

            /**
             * @brief Add a Dataset to the end of this Multiset.
             */
            void push_back(const Dataset2D& data);

            /**
             * @brief Add a Dataset to the end of this Multiset.
             */
            void push_back(const Dataset2D&& data);

            /**
             * @brief Impose a limit on the y-axis. All data lying outside this range will be removed.
             */
            void ylimits(const Limit& limit);

            /**
             * @brief Impose a limit on the y-axis. All data lying outside this range will be removed.
             */
            void ylimits(double min, double max);

            /**
             * @brief Save this Multiset at the given location.
             *        All constituent Datasets will be saved in a folder with the specified name. 
             */
            void save(const io::File& path) const;

            /**
             * @brief Read-only iterator.
             */
            const std::vector<Dataset2D>::const_iterator begin() const;

            /**
             * @brief Read-only iterator.
             */
            const std::vector<Dataset2D>::const_iterator end() const;

            /**
             * @brief Read-write iterator.
             */
            std::vector<Dataset2D>::iterator begin();

            /**
             * @brief Read-write iterator.
             */
            std::vector<Dataset2D>::iterator end();

            std::vector<Dataset2D> data;
            std::unordered_map<std::string, unsigned int> names;

        private:
            /**
             * @brief Read a saved Multiset.
             */
            void read(const io::Folder& path);
    };
}