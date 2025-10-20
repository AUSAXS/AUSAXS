// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/Definitions.h>

#include <unordered_map>
#include <functional>

namespace ausaxs::api {
    /**
     * @brief A simple storage for C++ objects that need to stay in memory for use in Python. 
     *        This avoids the need for C-style memory management despite being bound by a C API.
     */
    struct ObjectStorage {
        struct StoredObject {
            void* ptr;
            std::function<void(void*)> deleter;
        };
        template<typename T> static int register_object(T&& obj);

        /**
         * @brief Get an object by its ID. Manual type casting is required.
         */
        template<typename T> static T* get_object(int id);

        /**
         * @brief Deregister an object, deleting it and freeing its memory.
         */
        static void deregister_object(int id);

        static inline int current_id = 1;
        static inline std::unordered_map<int, StoredObject> storage;
    };

    template<typename T> 
    int ObjectStorage::register_object(T&& obj) {
        int id = current_id++;
        T* ptr = new T(std::move(obj));
        storage.emplace(id, StoredObject{
            .ptr=static_cast<void*>(ptr), 
            .deleter=[](void* p) { delete static_cast<T*>(p); }
        });
        return id;
    }

    template<typename T>
    inline T* ObjectStorage::get_object(int id) {
        auto it = storage.find(id);
        if (it != storage.end()) {
            return static_cast<T*>(it->second.ptr);
        }
        return nullptr;
    }

    inline void ObjectStorage::deregister_object(int id) {
        auto it = storage.find(id);
        if (it != storage.end()) {
            it->second.deleter(it->second.ptr);
            storage.erase(it);
        }
    }
}

extern "C" API void deallocate(int object_id, int* status);