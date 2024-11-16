#include <catch2/catch_test_macros.hpp>

#include <utility/StringUtils.h>

using namespace ausaxs;
using namespace utility;

TEST_CASE("remove_all") {
    SECTION("empty") {
        std::string result = utility::remove_all("", "");
        REQUIRE(result == "");
    }

    SECTION("single") {
        std::string result = utility::remove_all("a", "a");
        REQUIRE(result == "");
    }

    SECTION("multiple") {
        std::string result = utility::remove_all("aaa", "a");
        REQUIRE(result == "");
    }

    SECTION("multiple in middle") {
        std::string result = utility::remove_all("aabaabaa", "a");
        REQUIRE(result == "bb");
    }

    SECTION("multiple at start") {
        std::string result = utility::remove_all("aaabaa", "a");
        REQUIRE(result == "b");
    }

    SECTION("multiple at end") {
        std::string result = utility::remove_all("aaabaa", "a");
        REQUIRE(result == "b");
    }

    SECTION("multiple at start and end") {
        std::string result = utility::remove_all("aaaabaaa", "a");
        REQUIRE(result == "b");
    }

    SECTION("multiple in middle and at start and end") {
        std::string result = utility::remove_all("aaabaaabaaabaa", "a");
        REQUIRE(result == "bbb");
    }
}

TEST_CASE("StringUtils::remove_spaces") {
    SECTION("empty") {
        std::string result = utility::remove_spaces("");
        REQUIRE(result == "");
    }

    SECTION("single") {
        std::string result = utility::remove_spaces(" ");
        REQUIRE(result == "");
    }

    SECTION("multiple") {
        std::string result = utility::remove_spaces("   ");
        REQUIRE(result == "");
    }

    SECTION("single at start") {
        std::string result = utility::remove_spaces(" abc");
        REQUIRE(result == "abc");
    }

    SECTION("single at end") {
        std::string result = utility::remove_spaces("abc ");
        REQUIRE(result == "abc");
    }

    SECTION("multiple at start") {
        std::string result = utility::remove_spaces("   abc");
        REQUIRE(result == "abc");
    }

    SECTION("multiple at end") {
        std::string result = utility::remove_spaces("abc   ");
        REQUIRE(result == "abc");
    }

    SECTION("multiple at start and end") {
        std::string result = utility::remove_spaces("   abc   ");
        REQUIRE(result == "abc");
    }

    SECTION("multiple in middle") {
        std::string result = utility::remove_spaces("a   b   c");
        REQUIRE(result == "abc");
    }
}

TEST_CASE("StringUtils::remove_quotation_marks") {
    SECTION("empty") {
        std::string result = utility::remove_quotation_marks("");
        REQUIRE(result == "");
    }

    SECTION("single") {
        std::string result = utility::remove_quotation_marks("\"");
        REQUIRE(result == "\"");
    }

    SECTION("multiple") {
        std::string result = utility::remove_quotation_marks("\"\"\"");
        REQUIRE(result == "\"");
    }

    SECTION("basic") {
        std::string result = utility::remove_quotation_marks("\"abc\"");
        REQUIRE(result == "abc");
    }

    SECTION("single at start") {
        std::string result = utility::remove_quotation_marks("\"abc");
        REQUIRE(result == "\"abc");
    }

    SECTION("single at end") {
        std::string result = utility::remove_quotation_marks("abc\"");
        REQUIRE(result == "abc\"");
    }

    SECTION("multiple at start and end") {
        std::string result = utility::remove_quotation_marks("\"\"\"abc\"\"\"");
        REQUIRE(result == "\"\"abc\"\"");
    }
}

TEST_CASE("StringUtils::to_lowercase") {
    SECTION("empty") {
        std::string result = utility::to_lowercase("");
        REQUIRE(result == "");
    }

    SECTION("single") {
        std::string result = utility::to_lowercase("A");
        REQUIRE(result == "a");
    }

    SECTION("multiple") {
        std::string result = utility::to_lowercase("AbC");
        REQUIRE(result == "abc");
    }
}

TEST_CASE("StringUtils::split") {
    SECTION("single delimiter") {
        SECTION("empty") {
            std::vector<std::string> result = utility::split("", ',');
            REQUIRE(result.size() == 0);
        }

        SECTION("single") {
            std::vector<std::string> result = utility::split("a", ',');
            REQUIRE(result.size() == 1);
            REQUIRE(result[0] == "a");
        }

        SECTION("multiple") {
            std::vector<std::string> result = utility::split("a,b,c", ',');
            REQUIRE(result.size() == 3);
            REQUIRE(result[0] == "a");
            REQUIRE(result[1] == "b");
            REQUIRE(result[2] == "c");
        }
    }

    SECTION("multiple delimiters") {
        SECTION("empty") {
            std::vector<std::string> result = utility::split("", ",");
            REQUIRE(result.size() == 0);
        }

        SECTION("single") {
            std::vector<std::string> result = utility::split("a", ",");
            REQUIRE(result.size() == 1);
            REQUIRE(result[0] == "a");
        }

        SECTION("multiple") {
            std::vector<std::string> result = utility::split("a,b,c", ",");
            REQUIRE(result.size() == 3);
            REQUIRE(result[0] == "a");
            REQUIRE(result[1] == "b");
            REQUIRE(result[2] == "c");
        }

        SECTION("multiple with spaces") {
            std::vector<std::string> result = utility::split("a, b, c", ", ");
            REQUIRE(result.size() == 3);
            REQUIRE(result[0] == "a");
            REQUIRE(result[1] == "b");
            REQUIRE(result[2] == "c");
        }

        SECTION("multiple with lots of stuff") {
            std::vector<std::string> result = utility::split("  ,aaba  ,, ,,bda,  ced,,,  ,", ", ");
            REQUIRE(result.size() == 3);
            REQUIRE(result[0] == "aaba");
            REQUIRE(result[1] == "bda");
            REQUIRE(result[2] == "ced");
        }

        SECTION("actual example") {
            std::vector<std::string> result = utility::split("0.009813      	0.00667934    	0.00133647    	1\n\r", ", \t\n\r");
            REQUIRE(result.size() == 4);
            REQUIRE(result[0] == "0.009813");
            REQUIRE(result[1] == "0.00667934");
            REQUIRE(result[2] == "0.00133647");
            REQUIRE(result[3] == "1");
        }
    }
}

TEST_CASE("StringUtils::join") {
    SECTION("empty") {
        std::string result = utility::join({}, ",");
        REQUIRE(result == "");
    }

    SECTION("single") {
        std::string result = utility::join({"a"}, ",");
        REQUIRE(result == "a");
    }

    SECTION("multiple") {
        std::string result = utility::join({"a", "b", "c"}, ",");
        REQUIRE(result == "a,b,c");
    }

    SECTION("multiple with spaces") {
        std::string result = utility::join({"a", " b", " c"}, ",");
        REQUIRE(result == "a, b, c");
    }
}

TEST_CASE("StringUtils::remove_all") {
    SECTION("empty") {
        REQUIRE(utility::remove_all("", "a") == "");
    }

    SECTION("single") {
        REQUIRE(utility::remove_all("a", "a") == "");
    }

    SECTION("multiple") {
        REQUIRE(utility::remove_all("aaa", "a") == "");
    }

    SECTION("multiple with other") {
        REQUIRE(utility::remove_all("ababab", "a") == "bbb");
    }

    SECTION("multiple with other 2") {
        REQUIRE(utility::remove_all("ababab", "b") == "aaa");
    }

    SECTION("real issues") {
        REQUIRE(utility::remove_all("TITLE \r", " \r") == "TITLE");
    }
}