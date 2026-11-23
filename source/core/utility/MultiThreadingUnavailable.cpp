// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/MultiThreading.h>

#include <utility>

using namespace ausaxs;
using namespace ausaxs::utility::multi_threading;

namespace ausaxs::utility::multi_threading {
    class ThreadPoolImpl {};
}

ThreadPool::ThreadPool(int /*unused*/) : impl(std::make_unique<ThreadPoolImpl>()) {}
ThreadPool::~ThreadPool() = default;

void ThreadPool::detach_task(std::function<void()> task) {task();} // NOLINT

void ThreadPool::detach_blocks(int first, int last, std::function<void(int, int)> task) { // NOLINT
    if (first < last) {task(first, last);}
}

void Tasks::wait() const {wait_();}

Tasks ThreadPool::submit_blocks(int first, int last, std::function<void(int, int)> task) {
    detach_blocks(first, last, std::move(task));
    return Tasks{.wait_ = [] {}};
}

void ThreadPool::wait() {}
std::vector<std::thread::id> ThreadPool::get_thread_ids() const {return {};} // NOLINT

observer_ptr<ThreadPool> utility::multi_threading::get_global_pool() {
    static ThreadPool pool(1);
    return &pool;
}
