#pragma once

#include <map>
#include <string>

#include <utility/Exceptions.h>
#include <utility/Utility.h>

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
        SimpleMap(std::map<std::string, V> map) : data(map) {}

        /**
         * @brief Get a value from the storage. 
         */
        const V& get(std::string key) const {
            if (data.find(key) == data.end()) {
                throw except::map_error("Key " + key + " not found in map");
            }
            return data.at(key);
        }

        /**
         * @brief Get a value from the storage. 
         */
        V& get(std::string key) {
            return const_cast<V&>(std::as_const(*this).get(key));
        }

        /**
         * @brief Insert a key-value pair into the storage. 
         */
        void insert(std::string key, V val) {
            std::string k2 = utility::to_lowercase(key);
            data.emplace(k2, val);
        }

        /**
         * @brief Check if this map contains the given key. 
         */
        bool contains(std::string key) const {
            return data.find(utility::to_lowercase(key)) != data.end();
        }

        std::map<std::string, V> data;
    };
}