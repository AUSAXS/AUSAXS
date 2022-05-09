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
         * @brief Add a Dataset to the end of this Multiset.
         */
        void push_back(const Dataset& data);

        /**
         * @brief Add a Dataset to the end of this Multiset.
         */
        void push_back(const Dataset&& data);

    private:
        std::vector<Dataset> data;
        std::map<std::string, unsigned int> names;
};