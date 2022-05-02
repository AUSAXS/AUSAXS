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

    private:
        std::vector<Dataset> data;
        std::map<std::string, unsigned int> names;
};