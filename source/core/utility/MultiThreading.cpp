// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/MultiThreading.h>
#include <settings/GeneralSettings.h>

#include <BS_thread_pool.hpp>
#include <utility>

using namespace ausaxs;
using namespace ausaxs::utility::multi_threading;

class ThreadPool::Impl {
    public:
        explicit Impl(int threads) : pool(threads) {}
        BS::light_thread_pool pool;
};

ThreadPool::ThreadPool(int threads) : impl(std::make_unique<Impl>(threads)) {}
ThreadPool::~ThreadPool() = default;

void ThreadPool::detach_task(std::function<void()> task) {
    impl->pool.detach_task(std::move(task));
}

void ThreadPool::detach_blocks(int first, int last, std::function<void(int, int)> task) {
    impl->pool.detach_blocks(first, last, std::move(task));
}

void ThreadPool::detach_loop(int first, int last, std::function<void(int)> task) {
    impl->pool.detach_loop(first, last, std::move(task));
}

ThreadPool::Tasks::Tasks(std::function<void()> wait) : wait_(std::move(wait)) {}
void ThreadPool::Tasks::wait() const {wait_();}

ThreadPool::Tasks ThreadPool::submit_blocks(int first, int last, std::function<void(int, int)> task) {
    auto futures = std::make_shared<BS::multi_future<void>>(impl->pool.submit_blocks(first, last, std::move(task)));
    return Tasks([futures] {futures->wait();});
}

void ThreadPool::wait() {impl->pool.wait();}
std::vector<std::thread::id> ThreadPool::get_thread_ids() const {return impl->pool.get_thread_ids();}

observer_ptr<ThreadPool> utility::multi_threading::get_global_pool() {
    static ThreadPool pool(settings::general::threads);
    return &pool;
}
