// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>
#include <utility/Exceptions.h>

#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <typeindex>
#include <unordered_map>

namespace ausaxs::api {
    /**
     * @brief A simple storage for C++ objects that need to stay in memory for use in Python. 
     *        This avoids the need for C-style memory management despite being bound by a C API.
     */
    struct ObjectStorage {
        struct StoredObject {
            void* ptr = nullptr;
            std::function<void(void*)> deleter;
            std::type_index type = typeid(void); // the registered type; an object is only ever handed out as this type
        };
        template<typename T> static int register_object(T&& obj);
        template<typename T> static int register_object(std::unique_ptr<T> obj);

        /**
         * @brief Get an object by its ID.
         * @return The object, or nullptr if no object with this ID exists or it was registered as a different type.
         */
        template<typename T> static T* get_object(int id);

        /**
         * @brief Get an object by an ID the caller passed in.
         * @param kind What the object should be, in the caller's terms (e.g. "molecule"). It is only used in the error message.
         * @throws except::invalid_argument if no object with this ID exists, or if it was registered as a different type.
         */
        template<typename T> static T& get(int id, std::string_view kind);

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
        T* ptr = new T(std::forward<T>(obj)); // constructed outside the lock; it cannot be observed yet
        std::scoped_lock lock(mutex);
        int id = current_id++;
        storage.emplace(id, StoredObject{
            .ptr=static_cast<void*>(ptr), 
            .deleter=[](void* p) { delete static_cast<T*>(p); },
            .type=typeid(T)
        });
        return id;
    }

    template<typename T>
    int ObjectStorage::register_object(std::unique_ptr<T> obj) {
        T* ptr = obj.release(); // take ownership
        std::scoped_lock lock(mutex);
        int id = current_id++;
        storage.emplace(id, StoredObject{ ptr,
            [](void* p){ delete static_cast<T*>(p); },
            typeid(T)
        });
        return id;
    }

    template<typename T>
    inline T* ObjectStorage::get_object(int id) {
        std::scoped_lock lock(mutex);
        auto it = storage.find(id);
        if (it != storage.end() && it->second.type == typeid(T)) {
            return static_cast<T*>(it->second.ptr);
        }
        return nullptr;
    }

    template<typename T>
    inline T& ObjectStorage::get(int id, std::string_view kind) {
        std::scoped_lock lock(mutex);
        auto it = storage.find(id);
        if (it == storage.end()) {
            throw except::invalid_argument("No object with id \"" + std::to_string(id) + "\" exists (expected: " + std::string(kind) + "). It may have been deallocated.");
        }
        if (it->second.type != typeid(T)) {
            throw except::invalid_argument("Object \"" + std::to_string(id) + "\" has the wrong type (expected: " + std::string(kind) + ").");
        }
        return *static_cast<T*>(it->second.ptr);
    }

    inline void ObjectStorage::deregister_object(int id) {
        StoredObject obj;
        {
            std::scoped_lock lock(mutex);
            auto it = storage.find(id);
            if (it == storage.end()) {return;}
            obj = std::move(it->second);
            storage.erase(it);
        }
        obj.deleter(obj.ptr);
    }
}

extern "C" API void deallocate(int object_id, int* status);