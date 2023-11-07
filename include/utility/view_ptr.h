#pragma once

#include <memory>

/**
 * @brief A view-pointer is just a wrapper around a raw pointer to denote that ownership is not transferred.
 *        There are no guarantees on the lifetime of the object pointed to by a view-pointer. This is the user's responsibility.
 */
template<typename T>
class view_ptr {
    public:
        view_ptr(const view_ptr<T>&) = default;

        view_ptr(view_ptr<T>&& other) = default;

        template<typename U, typename = std::enable_if_t<std::is_base_of_v<T, U>>>
        view_ptr(U* ptr) : ptr(ptr) {}

        template<typename U, typename = std::enable_if_t<std::is_base_of_v<T, U>>>
        view_ptr(U& ref) : ptr(&ref) {}

        template<typename U, typename = std::enable_if_t<std::is_base_of_v<T, U>>>
        view_ptr(const std::unique_ptr<U>& ptr) : ptr(ptr.get()) {}

        template<typename U, typename = std::enable_if_t<std::is_base_of_v<T, U>>>
        view_ptr(const std::shared_ptr<U>& ptr) : ptr(ptr.get()) {}

        view_ptr& operator=(const view_ptr<T>& other) {
            ptr = other.ptr;
            return *this;
        }

        view_ptr& operator=(T* ptr) {
            this->ptr = ptr;
            return *this;
        }

        view_ptr& operator=(T& ref) {
            ptr = &ref;
            return *this;
        }

        T* get() {
            return ptr;
        }

        const T* get() const {
            return ptr;
        }

        T* operator->() const {
            return ptr;
        }

        T& operator*() const {
            return *ptr;
        }

        view_ptr<T> operator=(view_ptr<T> ptr) {
            this->ptr = ptr.ptr;
            return *this;
        }

        bool operator==(const view_ptr<T>& other) const {
            return ptr == other.ptr;
        }

    private:
        T* ptr;
};