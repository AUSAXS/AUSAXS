#pragma once

#include <map>
#include <string>
#include <iostream>

#include <utility/Exceptions.h>
#include <utility/Utility.h>

namespace saxs {
    namespace detail {
        /**
         * @brief A simple case-insensitive map. 
         */
        template<typename V>
        struct SimpleMap {
            /**
             * @brief Create a new empty SimpleMap.
             */
            SimpleMap() {}

            /**
             * @brief Create a new SimpleMap from a std::map.
             */
            SimpleMap(std::map<std::string, V> map) {
                for (auto& [key, value] : map) {
                    insert(key, value);
                }
            }

            /**
             * @brief Get a value from the storage. 
             */
            virtual const V& get(const std::string& key) const {
                std::string k2 = utility::to_lowercase(key);
                if (data.find(k2) == data.end()) {
                    throw except::map_error("Error in SimpleMap::get: Key " + k2 + " not found in map");
                }
                return data.at(k2);
            }

            /**
             * @brief Get a value from the storage. 
             */
            V& get(const std::string& key) {
                return const_cast<V&>(std::as_const(*this).get(key));
            }

            /**
             * @brief Insert a key-value pair into the storage. 
             */
            void insert(const std::string& key, V val) {
                std::string k2 = utility::to_lowercase(key);
                data.emplace(k2, val);
            }

            /**
             * @brief Check if this map contains the given key. 
             */
            bool contains(const std::string& key) const {
                return data.find(utility::to_lowercase(key)) != data.end();
            }

            typename std::map<std::string, V>::const_iterator begin() const {
                return data.begin();
            }

            typename std::map<std::string, V>::const_iterator end() const {
                return data.end();
            }

            typename std::map<std::string, V>::iterator begin() {
                return data.begin();
            }

            typename std::map<std::string, V>::iterator end() {
                return data.end();
            }

            std::map<std::string, V> data;
        };

        /**
         * @brief A simple extension of a SimpleMap, specializing it for residue storage. 
         *        The only difference is that all keys containing "h" are automatically mapped to 0. 
         */
        struct SimpleResidueMap : SimpleMap<unsigned int> {
            /**
             * @brief Create a new empty SimpleResidueMap.
             */
            SimpleResidueMap() : SimpleMap() {
                insert("h", 0);
            }

            /**
             * @brief Create a new SimpleResidueMap from a std::map.
             */
            SimpleResidueMap(std::map<std::string, unsigned int> map) : SimpleMap(map) {
                insert("h", 0);
            }

            /**
             * @brief Get a value from the storage. 
             */
            const unsigned int& get(const std::string& key) const override {
                std::string k2 = utility::to_lowercase(key);
                if (k2[0] == 'h') {
                    return data.at("h");
                }
                if (data.find(k2) == data.end()) {
                    throw except::map_error("Error in SimpleMap::get: Key " + k2 + " not found in map");
                }
                return data.at(k2);
            }
        };
    }
}