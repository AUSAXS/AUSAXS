// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <dataset/NamedWrapper.h>
#include <dataset/Dataset.h>
#include <dataset/Dataset2D.h>
#include <dataset/SimpleDataset.h>
#include <utility/Exceptions.h>
#include <io/File.h>

#include <string>
#include <fstream>

using namespace ausaxs;

template<typename T>
void NamedWrapper<T>::set_default_names() {
    names.resize(dataset.size_cols());
    for (unsigned int i = 0; i < dataset.size_cols(); i++) {
        names[i] = "col_" + std::to_string(i);
    }
}

template<typename T>
void NamedWrapper<T>::set_col_names(const std::vector<std::string>& new_names) {
    if (new_names.size() != dataset.size_cols()) {
        throw except::invalid_operation(
            "NamedWrapper::set_col_names: Number of names does not match number of columns. "
            "(" + std::to_string(new_names.size()) + " != " + std::to_string(dataset.size_cols()) + ")"
        );
    }
    names = new_names;
}

template<typename T>
void NamedWrapper<T>::set_col_names(unsigned int i, const std::string& name) {
    names[i] = name;
}

template<typename T>
std::vector<std::string> NamedWrapper<T>::get_col_names() const {
    return names;
}

template<typename T>
std::string NamedWrapper<T>::get_col_names(unsigned int i) const {
    return names[i];
}

template<typename T>
bool NamedWrapper<T>::is_named() const noexcept {
    for (unsigned int i = 0; i < dataset.size_cols(); i++) {
        if (names[i] != "col_" + std::to_string(i)) {
            return true;
        }
    }
    return false;
}

template<typename T>
MutableColumn<double> NamedWrapper<T>::col(std::string_view column) {
    for (unsigned int i = 0; i < names.size(); ++i) {
        if (names[i] == column) {
            return dataset.col(i);
        }
    }
    throw except::invalid_argument("NamedWrapper::col: Column \"" + std::string(column) + "\" not found.");
}

template<typename T>
const ConstColumn<double> NamedWrapper<T>::col(std::string_view column) const {
    for (unsigned int i = 0; i < names.size(); ++i) {
        if (names[i] == column) {
            return dataset.col(i);
        }
    }
    throw except::invalid_argument("NamedWrapper::col: Column \"" + std::string(column) + "\" not found.");
}

template<typename T>
NamedWrapper<Dataset> NamedWrapper<T>::select_columns(const std::vector<std::string>& cols) const {
    std::vector<unsigned int> indices(cols.size());
    std::transform(cols.begin(), cols.end(), indices.begin(), [this](const std::string& name) {
        for (unsigned int i = 0; i < names.size(); i++) {
            if (names[i] == name) {return i;}
        }
        throw except::invalid_argument("NamedWrapper::select_columns: Column \"" + name + "\" not found.");
    });
    
    Dataset new_dataset = dataset.select_columns(indices);
    std::vector<std::string> col_names(cols.size());
    for (unsigned int i = 0; i < cols.size(); i++) {
        col_names[i] = names[indices[i]];
    }
    return NamedWrapper<Dataset>(std::move(new_dataset), col_names);
}

template<typename T>
void NamedWrapper<T>::save(const io::File& path, const std::string& header) const {
    path.directory().create();

    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::ios_base::failure("NamedWrapper::save: Could not open file \"" + path.str() + "\"");
    }

    // write header
    if (!header.empty()) {
        output << header << std::endl;
    }

    // write column titles
    if (names.size() < dataset.size_cols()) {
        throw except::unexpected(
            "NamedWrapper::save: Number of column names (" + std::to_string(names.size()) + ") "
            "does not match number of columns (" + std::to_string(dataset.size_cols()) + ")."
        );
    }
    for (unsigned int j = 0; j < dataset.size_cols(); j++) {
        output << std::left << std::setw(16) << names[j] << "\t";
    }
    output << std::endl;

    // write data
    for (unsigned int i = 0; i < dataset.size_rows(); i++) {
        for (unsigned int j = 0; j < dataset.size_cols()-1; j++) {
            output << std::left << std::setw(16) << std::setprecision(8) << std::scientific << dataset.index(i, j) << "\t";
        }
        output << dataset.index(i, dataset.size_cols()-1) << "\n";
    }
    output.close();
}

template<typename T>
const T& NamedWrapper<T>::get_dataset() const noexcept {return dataset;}

template<typename T>
T& NamedWrapper<T>::get_dataset() noexcept {return dataset;}

template struct ausaxs::NamedWrapper<Dataset>;
template struct ausaxs::NamedWrapper<Dataset2D>;
template struct ausaxs::NamedWrapper<SimpleDataset>;