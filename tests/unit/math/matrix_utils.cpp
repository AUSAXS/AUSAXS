#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/MatrixUtils.h>
#include <math/Matrix.h>
#include <math/Vector3.h>

#include <numbers>

using namespace ausaxs;

TEST_CASE("matrix::identity") {
    SECTION("1x1") {
        Matrix<double> I = matrix::identity(1);
        REQUIRE(I.N == 1);
        REQUIRE(I.M == 1);
        REQUIRE(I[0][0] == 1);
    }

    SECTION("2x2") {
        Matrix<double> I = matrix::identity(2);
        REQUIRE(I == Matrix<double>{{1, 0}, {0, 1}});
    }

    SECTION("3x3") {
        Matrix<double> I = matrix::identity(3);
        REQUIRE(I == Matrix<double>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    }

    SECTION("4x4") {
        Matrix<double> I = matrix::identity(4);
        REQUIRE(I.N == 4);
        REQUIRE(I.M == 4);
        for (unsigned int i = 0; i < 4; ++i) {
            for (unsigned int j = 0; j < 4; ++j) {
                if (i == j) {
                    REQUIRE(I[i][j] == 1);
                } else {
                    REQUIRE(I[i][j] == 0);
                }
            }
        }
    }
}

TEST_CASE("matrix::rotation_matrix (Euler angles)") {
    SECTION("rotation around x-axis") {
        Matrix R = matrix::rotation_matrix(std::numbers::pi/2, 0., 0.);
        REQUIRE(R == Matrix{{1, 0, 0}, {0, 0, -1}, {0, 1, 0}});
    }

    SECTION("rotation around y-axis") {
        Matrix R = matrix::rotation_matrix(0., std::numbers::pi/2, 0.);
        REQUIRE(R == Matrix{{0, 0, 1}, {0, 1, 0}, {-1, 0, 0}});
    }

    SECTION("rotation around z-axis") {
        Matrix R = matrix::rotation_matrix(0., 0., std::numbers::pi/2);
        REQUIRE(R == Matrix{{0, -1, 0}, {1, 0, 0}, {0, 0, 1}});
    }

    SECTION("no rotation") {
        Matrix R = matrix::rotation_matrix(0., 0., 0.);
        REQUIRE(R == matrix::identity(3));
    }

    SECTION("rotation from Vector3") {
        Vector3<double> angles = {std::numbers::pi/2, 0., 0.};
        Matrix R = matrix::rotation_matrix(angles);
        REQUIRE(R == Matrix{{1, 0, 0}, {0, 0, -1}, {0, 1, 0}});
    }
}

TEST_CASE("matrix::rotation_matrix (axis-angle)") {
    SECTION("rotation around x-axis") {
        Vector3<double> axis = {1, 0, 0};
        Matrix R = matrix::rotation_matrix(axis, std::numbers::pi/2);
        
        Vector3<double> v = {0, 1, 0};
        Vector3<double> result = R * v;
        REQUIRE(result == Vector3<double>{0, 0, 1});
    }

    SECTION("rotation around y-axis") {
        Vector3<double> axis = {0, 1, 0};
        Matrix R = matrix::rotation_matrix(axis, std::numbers::pi/2);
        
        Vector3<double> v = {1, 0, 0};
        Vector3<double> result = R * v;
        REQUIRE(result == Vector3<double>{0, 0, -1});
    }

    SECTION("rotation around z-axis") {
        Vector3<double> axis = {0, 0, 1};
        Matrix R = matrix::rotation_matrix(axis, std::numbers::pi/2);
        
        Vector3<double> v = {1, 0, 0};
        Vector3<double> result = R * v;
        REQUIRE(result == Vector3<double>{0, 1, 0});
    }

    SECTION("no rotation") {
        Vector3<double> axis = {1, 0, 0};
        Matrix R = matrix::rotation_matrix(axis, 0.0);
        REQUIRE(R == matrix::identity(3));
    }
}

TEST_CASE("matrix::rotation_matrix orthonormality") {
    SECTION("Euler angles") {
        for (int i = 0; i < 10; i++) {
            double alpha = (rand() % 100) * std::numbers::pi / 100;
            double beta = (rand() % 100) * std::numbers::pi / 100;
            double gamma = (rand() % 100) * std::numbers::pi / 100;
            
            Matrix R = matrix::rotation_matrix(alpha, beta, gamma);
            Matrix Rt = R.T();
            Matrix I = R * Rt;
            
            REQUIRE(I == matrix::identity(3));
        }
    }

    SECTION("axis-angle") {
        for (int i = 0; i < 10; i++) {
            Vector3<double> axis = {
                static_cast<double>(rand() % 100),
                static_cast<double>(rand() % 100),
                static_cast<double>(rand() % 100)
            };
            double angle = (rand() % 100) * std::numbers::pi / 100;
            
            Matrix R = matrix::rotation_matrix(axis, angle);
            Matrix Rt = R.T();
            Matrix I = R * Rt;
            
            REQUIRE(I == matrix::identity(3));
        }
    }
}

TEST_CASE("vector3::generate_basis") {
    SECTION("simple case") {
        Vector3<double> v = {2, 0, 0};
        auto [v1, v2, v3] = vector3::generate_basis(v);
        
        REQUIRE(v1 == Vector3<double>{1, 0, 0});
        REQUIRE_THAT(v1.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(v2.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(v3.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(v1.dot(v2), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(v2.dot(v3), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(v3.dot(v1), Catch::Matchers::WithinAbs(0.0, 1e-10));
    }

    SECTION("arbitrary vector") {
        Vector3<double> v = {1, 2, 3};
        auto [v1, v2, v3] = vector3::generate_basis(v);
        
        REQUIRE_THAT(v1.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(v2.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(v3.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(v1.dot(v2), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(v2.dot(v3), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(v3.dot(v1), Catch::Matchers::WithinAbs(0.0, 1e-10));
    }
}
