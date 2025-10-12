#include <catch2/catch_test_macros.hpp>

#include <container/ThreadLocalWrapper.h>
#include <container/Container1D.h>
#include <settings/GeneralSettings.h>
#include <utility/MultiThreading.h>

using namespace ausaxs;
using namespace container;

TEST_CASE("ThreadLocalWrapper::ThreadLocalWrapper") {
    SECTION("default") {
        ThreadLocalWrapper<int> wrapper;
        CHECK(wrapper.get() == 0);
        CHECK(wrapper.size() == settings::general::threads+1);
    }

    SECTION("Args&&...") {
        SECTION("int") {
            ThreadLocalWrapper<int> wrapper(5);
            CHECK(wrapper.get() == 5);
            CHECK(wrapper.size() == settings::general::threads+1);
        }

        SECTION("container1D") {
            ThreadLocalWrapper<Container1D<double>> wrapper(std::vector<double>{0.5, 1, 1.5, 2});
            CHECK(wrapper.size() == settings::general::threads+1);
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
    for (unsigned int i = 0; i < 100; ++i) {
        pool->detach_task([&wrapper](){wrapper.get() += 1;});
    }
    pool->wait();
    auto merged = wrapper.merge();
    CHECK(merged == 100);
}