// #include <catch2/catch_test_macros.hpp>
// #include <catch2/matchers/catch_matchers_floating_point.hpp>

// #include <crystal/CrystalScattering.h>
// #include <crystal/miller/AllMillers.h>
// #include <crystal/miller/FibonacciMillers.h>

// using namespace crystal;

// TEST_CASE("millers_generation") {
    // std::shared_ptr<MillerGenerationStrategy> miller_strategy = std::make_shared<AllMillers>(10, 10, 10);
    // std::shared_ptr<MillerGenerationStrategy> miller_strategy = std::make_shared<FibonacciMillers>(10, 10, 10);
    // std::vector<Miller> millers = miller_strategy->generate();

    // crystal::FibonacciMillers fib(10, 10, 10);
    // std::vector<crystal::Miller> millers = fib.generate();
    // fib.generate_fibonacci_sphere(100);
// }

// TEST_CASE("debug") {
//     setting::crystal::mgc = setting::crystal::Fibonacci;
//     Grid grid = Grid({10, -10, 10, -10, 10, -10, 1}, 1);
//     CrystalScattering cs(grid);
//     // cs.calculate();
// }