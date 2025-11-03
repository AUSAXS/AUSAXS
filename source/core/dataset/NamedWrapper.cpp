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
    names.resize(this->size_cols());
    for (unsigned int i = 0; i < this->size_cols(); i++) {
        names[i] = "col_" + std::to_string(i);
    }
}

template<typename T>
void NamedWrapper<T>::set_col_names(const std::vector<std::string>& new_names) {
    if (new_names.size() != this->size_cols()) {
        throw except::invalid_operation(
            "NamedWrapper::set_col_names: Number of names does not match number of columns. "
            "(" + std::to_string(new_names.size()) + " != " + std::to_string(this->size_cols()) + ")"
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
    for (unsigned int i = 0; i < this->size_cols(); i++) {
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
            return this->col(i);
        }
    }
    throw except::invalid_argument("NamedWrapper::col: Column \"" + std::string(column) + "\" not found.");
}

template<typename T>
const ConstColumn<double> NamedWrapper<T>::col(std::string_view column) const {
    for (unsigned int i = 0; i < names.size(); ++i) {
        if (names[i] == column) {
            return this->col(i);
        }
    }
    throw except::invalid_argument("NamedWrapper::col: Column \"" + std::string(column) + "\" not found.");
}

template<typename T>
NamedWrapper<Dataset> NamedWrapper<T>::select_columns(std::initializer_list<std::string_view> cols) const {
    std::vector<unsigned int> indices(cols.size());
    std::transform(cols.begin(), cols.end(), indices.begin(), [this](std::string_view name) {
        for (unsigned int i = 0; i < names.size(); i++) {
            if (names[i] == name) {return i;}
        }
        throw except::invalid_argument("NamedWrapper::select_columns: Column \"" + std::string(name) + "\" not found.");
    });
    
    NamedWrapper<Dataset> new_dataset = this->select_columns(indices);
    std::vector<std::string> col_names(cols.size());
    for (unsigned int i = 0; i < cols.size(); i++) {
        col_names[i] = names[indices[i]];
    }
    new_dataset.names = std::move(col_names);
    return new_dataset;
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
    if (names.size() < this->size_cols()) {
        throw except::unexpected(
            "NamedWrapper::save: Number of column names (" + std::to_string(names.size()) + ") "
            "does not match number of columns (" + std::to_string(this->size_cols()) + ")."
        );
    }
    for (unsigned int j = 0; j < this->size_cols(); j++) {
        output << std::left << std::setw(16) << names[j] << "\t";
    }
    output << std::endl;

    // write data
    for (unsigned int i = 0; i < this->size_rows(); i++) {
        for (unsigned int j = 0; j < this->size_cols()-1; j++) {
            output << std::left << std::setw(16) << std::setprecision(8) << std::scientific << this->index(i, j) << "\t";
        }
        output << this->index(i, this->size_cols()-1) << "\n";
    }
    output.close();
}

template struct ausaxs::NamedWrapper<Dataset>;
template struct ausaxs::NamedWrapper<Dataset2D>;
template struct ausaxs::NamedWrapper<SimpleDataset>;