#include <catch2/catch_test_macros.hpp>

#include <dataset/SimpleDataset.h>

TEST_CASE("SimpleDataset::SimpleDataset") {
    SECTION("default constructor") {}
    SECTION("Dataset&") {}
    SECTION("unsigned int") {}
    SECTION("std::vector<double>, std::vector<double>, std::vector<double>") {}
    SECTION("std::vector<double>, std::vector<double>") {}
    SECTION("std::vector<double>, std::vector<double>, std::string, std::string") {}
    SECTION("io::ExistingFile") {}
}

TEST_CASE("SimpleDataset::yerr") {}
TEST_CASE("SimpleDataset::load") {}
TEST_CASE("SimpleDataset::reduce") {}
TEST_CASE("SimpleDataset::operator=") {}
TEST_CASE("SimpleDataset::operator==") {}
TEST_CASE("SimpleDataset::span_x") {}
TEST_CASE("SimpleDataset::span_y") {}
TEST_CASE("SimpleDataset::get_xlimits") {}
TEST_CASE("SimpleDataset::get_ylimits") {}
TEST_CASE("SimpleDataset::span_y_positive") {}
TEST_CASE("SimpleDataset::push_back") {}
TEST_CASE("SimpleDataset::normalize") {}
TEST_CASE("SimpleDataset::scale_errors") {}
TEST_CASE("SimpleDataset::scale_y") {}
TEST_CASE("SimpleDataset::simulate_noise") {}
TEST_CASE("SimpleDataset::simulate_errors") {}
TEST_CASE("SimpleDataset::get_point") {}
TEST_CASE("SimpleDataset::find_minimum") {}
TEST_CASE("SimpleDataset::rebin") {}
TEST_CASE("SimpleDataset::generate_random_data") {}
TEST_CASE("SimpleDataset::mean") {}
TEST_CASE("SimpleDataset::weighted_mean") {}
TEST_CASE("SimpleDataset::std") {}
TEST_CASE("SimpleDataset::weighted_mean_error") {}
TEST_CASE("SimpleDataset::remove_consecutive_duplicates") {}