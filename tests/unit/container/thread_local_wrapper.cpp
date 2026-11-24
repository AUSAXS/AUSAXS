#include <catch2/catch_test_macros.hpp>

#include <container/Container1D.h>
#include <container/ThreadLocalWrapper.h>
#include <settings/GeneralSettings.h>
#include <utility/MultiThreading.h>

#include <atomic>

using namespace ausaxs;
using namespace container;

TEST_CASE("ThreadLocalWrapper::ThreadLocalWrapper") {
    SECTION("default") {
        ThreadLocalWrapper<int> wrapper;
        CHECK(wrapper.get() == 0);
        CHECK(static_cast<int>(wrapper.size()) == settings::general::threads+1);
    }

    SECTION("Args&&...") {
        SECTION("int") {
            ThreadLocalWrapper<int> wrapper(5);
            CHECK(wrapper.get() == 5);
            CHECK(static_cast<int>(wrapper.size()) == settings::general::threads+1);
        }

        SECTION("container1D") {
            ThreadLocalWrapper<Container1D<double>> wrapper(std::vector<double>{0.5, 1, 1.5, 2});
            CHECK(static_cast<int>(wrapper.size()) == settings::general::threads+1);
            for (const auto& t : wrapper.get_all()) {
                CHECK(t.get().index(0) == 0.5);
                CHECK(t.get().index(1) == 1);
                CHECK(t.get().index(2) == 1.5);
                CHECK(t.get().index(3) == 2);
            }
        }
    }
}

TEST_CASE("ThreadLocalWrapper::ops") {
    const auto& pool = utility::multi_threading::get_global_pool();
    ThreadLocalWrapper<int> wrapper(0);
    for (int i = 0; i < 100; ++i) {
        pool->detach_task([&wrapper](){wrapper.get() += 1;});
    }
    pool->wait();
    auto merged = wrapper.merge();
    CHECK(merged == 100);
}

TEST_CASE("ThreadPool task operations") {
    utility::multi_threading::ThreadPool pool(2);
    std::atomic<int> sum{0};
    pool.detach_task([&] {sum += 1;});
    pool.detach_loop(0, 4, [&](int i) {sum += i;});
    pool.detach_blocks(3, 5, [&](int first, int last) {
        for (int i = first; i < last; ++i) {sum += i;}
    });
    pool.wait();
    CHECK(sum == 14);

    std::atomic<int> block_sum{0};
    pool.submit_blocks(0, 4, [&](int first, int last) {
        for (int i = first; i < last; ++i) {block_sum += i;}
    }).wait();
    CHECK(block_sum == 6);
    CHECK(pool.get_thread_ids().size() == 2);
}