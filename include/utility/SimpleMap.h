#pragma once

#include <utility/StringUtils.h>
#include <utility/Exceptions.h>

#include <unordered_map>
#include <string>

namespace saxs {
    namespace detail {
        /**
         * @brief A simple case-insensitive hashmap. 
         */
        template<typename V>
        struct SimpleMap {
            /**
             * @brief Create a new empty SimpleMap.
             */
            SimpleMap() = default;

            /**
             * @brief Create a new SimpleMap from a std::map.
             */
            SimpleMap(std::unordered_map<std::string, V> map) {
                for (auto& [key, value] : map) {
                    insert(key, value);
                }
            }

            /**
             * @brief Get a value from the storage. 
             */
            virtual V get(const std::string& key) const {
                std::string k2 = utility::to_lowercase(key);
                if (!data.contains(k2)) {
                    throw except::map_error("SimpleMap::get: Key " + k2 + " not found in map");
                }
                return data.at(k2);
            }

            /**
             * @brief Insert a key-value pair into the storage. 
             */
            void insert(const std::string& key, V val) {
                std::string k2 = utility::to_lowercase(key);
                data.emplace(k2, val);
            }

            const std::unordered_map<std::string, V>& get_map() const {
                return data;
            }

            /**
             * @brief Check if this map contains the given key. 
             */
            bool contains(const std::string& key) const {
                return data.contains(utility::to_lowercase(key));
            }

            typename std::unordered_map<std::string, V>::const_iterator begin() const {
                return data.begin();
            }

            typename std::unordered_map<std::string, V>::const_iterator end() const {
                return data.end();
            }

            typename std::unordered_map<std::string, V>::iterator begin() {
                return data.begin();
            }

            typename std::unordered_map<std::string, V>::iterator end() {
                return data.end();
            }

            std::unordered_map<std::string, V> data;
        };
    }
}