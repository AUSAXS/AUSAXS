// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

#include <unordered_map>
#include <functional>
#include <memory>
#include <mutex>

namespace ausaxs::api {
    /**
     * @brief A simple storage for C++ objects that need to stay in memory for use in Python. 
     *        This avoids the need for C-style memory management despite being bound by a C API.
     */
    struct ObjectStorage {
        struct StoredObject {
            void* ptr = nullptr;
            std::function<void(void*)> deleter;
        };
        template<typename T> static int register_object(T&& obj);
        template<typename T> static int register_object(std::unique_ptr<T> obj);

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
        static inline std::mutex mutex;
    };

    template<typename T> 
    int ObjectStorage::register_object(T&& obj) {
        T* ptr = new T(std::move(obj)); // constructed outside the lock; it cannot be observed yet
        std::lock_guard lock(mutex);
        int id = current_id++;
        storage.emplace(id, StoredObject{
            .ptr=static_cast<void*>(ptr), 
            .deleter=[](void* p) { delete static_cast<T*>(p); }
        });
        return id;
    }

    template<typename T>
    int ObjectStorage::register_object(std::unique_ptr<T> obj) {
        T* ptr = obj.release(); // take ownership
        std::lock_guard lock(mutex);
        int id = current_id++;
        storage.emplace(id, StoredObject{ ptr,
            [](void* p){ delete static_cast<T*>(p); }});
        return id;
    }

    template<typename T>
    inline T* ObjectStorage::get_object(int id) {
        std::lock_guard lock(mutex);
        auto it = storage.find(id);
        if (it != storage.end()) {
            return static_cast<T*>(it->second.ptr);
        }
        return nullptr;
    }

    inline void ObjectStorage::deregister_object(int id) {
        StoredObject obj;
        {
            std::lock_guard lock(mutex);
            auto it = storage.find(id);
            if (it == storage.end()) {return;}
            obj = std::move(it->second);
            storage.erase(it);
        }
        obj.deleter(obj.ptr);
    }
}

extern "C" API void deallocate(int object_id, int* status);