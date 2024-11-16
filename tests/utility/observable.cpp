#include <catch2/catch_test_macros.hpp>
#include <utility/Observable.h>

#include <iostream>

using namespace ausaxs;

TEST_CASE("Observables & Observers") {
    utility::Observable<int> observable;

    double observed_value = 0;
    {
        auto observer = observable.make_observer();

        observer->on_notify = [&observed_value] (const int& value) {
            observed_value = value;
        };

        observable.notify(1);
        REQUIRE(observed_value == 1);

        observable.notify(2);
        REQUIRE(observed_value == 2);

        observable.notify(3);
        REQUIRE(observed_value == 3);
    }

    observable.notify(4);
    REQUIRE(observed_value == 3);

    observable.notify(5);
    REQUIRE(observed_value == 3);

    observable.notify(6);
    REQUIRE(observed_value == 3);
}