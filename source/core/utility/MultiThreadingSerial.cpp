// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/MultiThreading.h>

#include <utility>

using namespace ausaxs;
using namespace ausaxs::utility::multi_threading;

class ThreadPool::Impl {};

ThreadPool::ThreadPool(int) : impl(std::make_unique<Impl>()) {}
ThreadPool::~ThreadPool() = default;

void ThreadPool::detach_task(std::function<void()> task) {task();}

void ThreadPool::detach_blocks(int first, int last, std::function<void(int, int)> task) {
    if (first < last) {task(first, last);}
}

void ThreadPool::detach_loop(int first, int last, std::function<void(int)> task) {
    for (int i = first; i < last; ++i) {task(i);}
}

ThreadPool::Tasks::Tasks(std::function<void()> wait) : wait_(std::move(wait)) {}
void ThreadPool::Tasks::wait() const {wait_();}

ThreadPool::Tasks ThreadPool::submit_blocks(int first, int last, std::function<void(int, int)> task) {
    detach_blocks(first, last, std::move(task));
    return Tasks([] {});
}

void ThreadPool::wait() {}
std::vector<std::thread::id> ThreadPool::get_thread_ids() const {return {};}

observer_ptr<ThreadPool> utility::multi_threading::get_global_pool() {
    static ThreadPool pool(1);
    return &pool;
}
