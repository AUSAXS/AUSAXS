// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <unordered_map>
#include <vector>

namespace ausaxs::api {
    struct ObjectStorage {
        static int register_object(void* obj);
        static int register_object(std::vector<void*> obj);
        template<typename T> static int register_object(T&& obj);
        template<typename T> static T* get_object(int id);
        static void unregister_object(int id);
        static inline std::unordered_map<int, void*> storage;
    };

    inline int ObjectStorage::register_object(void* obj) {
        static int current_id = 0;
        int id = current_id++;
        storage[id] = obj;
        return id;
    }

    template<typename T> 
    int ObjectStorage::register_object(T&& obj) {
        return register_object(static_cast<void*>(new T(std::forward<T>(obj))));
    }

    template<typename T>
    inline T* ObjectStorage::get_object(int id) {
        auto it = storage.find(id);
        if (it != storage.end()) {
            return static_cast<T*>(it->second);
        }
        return nullptr;
    }

    inline void ObjectStorage::unregister_object(int id) {
        storage.erase(id);
    }
}