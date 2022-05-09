#pragma once

#include <vector>
#include <string>

#include <utility/Dataset.h>

class Multiset {
    public:
        Multiset() {}

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

    private:
        std::vector<Dataset> data;
        std::map<std::string, unsigned int> names;
};