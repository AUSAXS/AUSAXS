#pragma once

#include <vector>
#include <string>

#include <utility/Dataset.h>

class Multiset {
    public:
        Multiset() {}

        Multiset(std::string path) {read(path);}

        explicit Multiset(unsigned int size) : data(size) {}

        explicit Multiset(const std::vector<Dataset>& data);

        explicit Multiset(const Dataset& data);

        Multiset(const Dataset& data1, const Dataset& data2);

        const Dataset& operator[](unsigned int i) const;
        Dataset& operator[](unsigned int i);

        const Dataset& get_data(std::string name) const;
        Dataset& get_data(std::string name);

        const Dataset& get_data(unsigned int i) const;
        Dataset& get_data(unsigned int i);

        /**
         * @brief Get the number of Datasets contained in this Multiset. 
         */
        size_t size() const;

        /**
         * @brief Check if this Multifram is empty.
         */
        bool empty() const;

        /**
         * @brief Add a Dataset to the end of this Multiset.
         */
        void push_back(const Dataset& data);

        /**
         * @brief Add a Dataset to the end of this Multiset.
         */
        void push_back(const Dataset&& data);

        /**
         * @brief Impose a limit on the y-axis. All data lying outside this range will be removed.
         */
        void ylimits(const Limit& limit) noexcept;

        /**
         * @brief Impose a limit on the y-axis. All data lying outside this range will be removed.
         */
        void ylimits(double min, double max) noexcept;

        /**
         * @brief Save this Multiset at the given location.
         *        All constituent Datasets will be saved in a folder with the specified name. 
         */
        void save(std::string path) const;

        /**
         * @brief Read-only iterator.
         */
		const std::vector<Dataset>::const_iterator begin() const;

        /**
         * @brief Read-only iterator.
         */
        const std::vector<Dataset>::const_iterator end() const;

        /**
         * @brief Read-write iterator.
         */
        std::vector<Dataset>::iterator begin();

        /**
         * @brief Read-write iterator.
         */
        std::vector<Dataset>::iterator end();

        std::vector<Dataset> data;
        std::map<std::string, unsigned int> names;

    private:
        /**
         * @brief Read a saved Multiset.
         */
        void read(std::string path);
};