#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hydrate/Hydration.h>
#include <hydrate/ExplicitHydration.h>
#include <hydrate/ImplicitHydration.h>
#include <hydrate/EmptyHydration.h>
#include <data/atoms/Water.h>

using namespace ausaxs;
using namespace ausaxs::hydrate;
using namespace ausaxs::data;

TEST_CASE("ExplicitHydration::ExplicitHydration") {
    SECTION("default constructor") {
        ExplicitHydration hydration;
        CHECK(hydration.waters.empty());
    }

    SECTION("constructor with vector copy") {
        std::vector<Water> waters = {
            Water({1, 2, 3}),
            Water({4, 5, 6}),
            Water({7, 8, 9})
        };
        ExplicitHydration hydration(waters);
        REQUIRE(hydration.waters.size() == 3);
        CHECK(hydration.waters[0] == Water({1, 2, 3}));
        CHECK(hydration.waters[1] == Water({4, 5, 6}));
        CHECK(hydration.waters[2] == Water({7, 8, 9}));
    }

    SECTION("constructor with vector move") {
        std::vector<Water> waters = {
            Water({1, 2, 3}),
            Water({4, 5, 6})
        };
        ExplicitHydration hydration(std::move(waters));
        REQUIRE(hydration.waters.size() == 2);
        CHECK(hydration.waters[0] == Water({1, 2, 3}));
        CHECK(hydration.waters[1] == Water({4, 5, 6}));
    }
}

TEST_CASE("ExplicitHydration::clear") {
    std::vector<Water> waters = {
        Water({1, 2, 3}),
        Water({4, 5, 6}),
        Water({7, 8, 9})
    };
    ExplicitHydration hydration(waters);
    REQUIRE(hydration.waters.size() == 3);

    hydration.clear();
    CHECK(hydration.waters.empty());
}

TEST_CASE("ExplicitHydration::clone") {
    std::vector<Water> waters = {
        Water({1, 2, 3}),
        Water({4, 5, 6})
    };
    ExplicitHydration hydration(waters);

    auto cloned = hydration.clone();
    REQUIRE(cloned != nullptr);

    auto* explicit_clone = dynamic_cast<ExplicitHydration*>(cloned.get());
    REQUIRE(explicit_clone != nullptr);
    REQUIRE(explicit_clone->waters.size() == 2);
    CHECK(explicit_clone->waters[0] == Water({1, 2, 3}));
    CHECK(explicit_clone->waters[1] == Water({4, 5, 6}));
}

TEST_CASE("ImplicitHydration::ImplicitHydration") {
    SECTION("constructor") {
        ImplicitHydration hydration;
        // Just verify it constructs successfully
        CHECK(true);
    }
}

TEST_CASE("ImplicitHydration::clear") {
    ImplicitHydration hydration;
    // clear() should throw for ImplicitHydration
    CHECK_THROWS_AS(hydration.clear(), std::runtime_error);
}

TEST_CASE("ImplicitHydration::clone") {
    ImplicitHydration hydration;
    // clone() should throw for ImplicitHydration
    CHECK_THROWS_AS(hydration.clone(), std::runtime_error);
}

TEST_CASE("EmptyHydration::EmptyHydration") {
    SECTION("constructor") {
        EmptyHydration hydration;
        // Just verify it constructs successfully
        CHECK(true);
    }
}

TEST_CASE("EmptyHydration::clear") {
    EmptyHydration hydration;
    // clear() should not throw and does nothing
    CHECK_NOTHROW(hydration.clear());
}

TEST_CASE("EmptyHydration::clone") {
    EmptyHydration hydration;
    auto cloned = hydration.clone();
    REQUIRE(cloned != nullptr);

    auto* empty_clone = dynamic_cast<EmptyHydration*>(cloned.get());
    CHECK(empty_clone != nullptr);
}

TEST_CASE("Hydration::create") {
    SECTION("create with empty constructor") {
        auto hydration = Hydration::create();
        REQUIRE(hydration != nullptr);

        auto* empty_hydration = dynamic_cast<EmptyHydration*>(hydration.get());
        CHECK(empty_hydration != nullptr);
    }

    SECTION("create with vector lvalue") {
        std::vector<Water> waters = {
            Water({1, 2, 3}),
            Water({4, 5, 6})
        };
        auto hydration = Hydration::create(waters);
        REQUIRE(hydration != nullptr);

        auto* explicit_hydration = dynamic_cast<ExplicitHydration*>(hydration.get());
        REQUIRE(explicit_hydration != nullptr);
        REQUIRE(explicit_hydration->waters.size() == 2);
        CHECK(explicit_hydration->waters[0] == Water({1, 2, 3}));
        CHECK(explicit_hydration->waters[1] == Water({4, 5, 6}));
    }

    SECTION("create with vector rvalue") {
        auto hydration = Hydration::create(std::vector<Water>{
            Water({7, 8, 9}),
            Water({10, 11, 12})
        });
        REQUIRE(hydration != nullptr);

        auto* explicit_hydration = dynamic_cast<ExplicitHydration*>(hydration.get());
        REQUIRE(explicit_hydration != nullptr);
        REQUIRE(explicit_hydration->waters.size() == 2);
        CHECK(explicit_hydration->waters[0] == Water({7, 8, 9}));
        CHECK(explicit_hydration->waters[1] == Water({10, 11, 12}));
    }
}
