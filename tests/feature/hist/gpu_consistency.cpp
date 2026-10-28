#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/state/Signaller.h>  // IWYU pragma: keep
#include <data/symmetry/CyclicSymmetry.h>
#include <gpu/GPULoader.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <rigidbody/BodySplitter.h>
#include <settings/All.h>

#include <hist/hist_test_helper.h>

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

/**
 * Every calculation that can be sent to the GPU is compared against the same calculation on the CPU.
 *
 * Without a GPU backend, as on the CI runners, the GPU kernel hands everything to the CPU calculator from its constructor,
 * so this then checks that fallback. With a backend and a device, it checks the device.
 */

using namespace ausaxs;
using namespace ausaxs::data;

namespace {
    /**
     * @brief Run @a f with the calculations sent to the GPU, if @a gpu.
     */
    template<typename F>
    auto on(bool gpu, F&& f) {
        settings::general::gpu = gpu;
        auto result = f();
        settings::general::gpu = false;
        return result;
    }

    template<typename T>
    std::vector<double> flat(const T& distribution) {return std::vector<double>(distribution.begin(), distribution.end());}

    /**
     * @brief The first moment of each bin, which carries the weighted bin centres without the noise of the sparsely populated ones.
     */
    std::vector<double> moments(const hist::ICompositeDistanceHistogram& h) {
        std::vector<double> m = h.get_weighted_counts();
        for (std::size_t i = 0; i < m.size(); ++i) {m[i] *= h.get_d_axis()[i];}
        return m;
    }

    // a pair right on a bin edge may round into the neighbouring bin on the device, whose sqrt differs from the CPU's in the last bits
    void check(const hist::ICompositeDistanceHistogram& cpu, const hist::ICompositeDistanceHistogram& gpu) {
        CHECK(compare_hist_approx(cpu.get_weighted_counts(), gpu.get_weighted_counts()));
        CHECK(compare_hist_approx(moments(cpu), moments(gpu)));

        const auto* cpu_exv = dynamic_cast<const hist::ICompositeDistanceHistogramExv*>(&cpu);
        const auto* gpu_exv = dynamic_cast<const hist::ICompositeDistanceHistogramExv*>(&gpu);
        REQUIRE((cpu_exv == nullptr) == (gpu_exv == nullptr));
        if (cpu_exv == nullptr) {
            CHECK(compare_hist_approx(flat(cpu.get_aa_counts()), flat(gpu.get_aa_counts())));
            CHECK(compare_hist_approx(flat(cpu.get_aw_counts()), flat(gpu.get_aw_counts())));
            CHECK(compare_hist_approx(flat(cpu.get_ww_counts()), flat(gpu.get_ww_counts())));
        } else {
            CHECK(compare_hist_approx(flat(cpu_exv->get_raw_aa_counts_by_ff()), flat(gpu_exv->get_raw_aa_counts_by_ff())));
            CHECK(compare_hist_approx(flat(cpu_exv->get_raw_aw_counts_by_ff()), flat(gpu_exv->get_raw_aw_counts_by_ff())));
            CHECK(compare_hist_approx(flat(cpu_exv->get_raw_ww_counts_by_ff()), flat(gpu_exv->get_raw_ww_counts_by_ff())));
            CHECK(compare_hist_approx(cpu_exv->get_total_raw_counts(), gpu_exv->get_total_raw_counts()));
        }
    }

    /**
     * @brief Compare a fresh @a Manager on the CPU and on the GPU.
     */
    template<typename Manager>
    void check_manager(const Molecule& protein) {
        auto cpu = on(false, [&] () {return Manager(&protein).calculate_all();});
        auto gpu = on(true,  [&] () {return Manager(&protein).calculate_all();});
        check(*cpu, *gpu);
    }

    /**
     * @brief 2epe split into three bodies, hydrated with @a waters, with a doubly repeated translational symmetry on the first.
     *        The hydration is passed in since generating it is random, and the CPU and GPU molecules must be identical.
     */
    Molecule symmetric_2epe(const std::vector<Water>& waters) {
        auto split = rigidbody::BodySplitter::split("tests/files/2epe.pdb", {40, 80});
        std::vector<Body> bodies;
        bodies.reserve(split.size_body());
        for (int i = 0; i < split.size_body(); ++i) {bodies.push_back(split.get_body(i));}
        bodies[0] = Body(std::as_const(bodies[0]).get_atoms(), waters);

        Molecule protein(std::move(bodies));
        protein.get_body(0).symmetry().add(std::make_unique<symmetry::CyclicSymmetry>(
            Vector3<double>{0, 0, 0}, Vector3<double>{-25, 0, 0}, Vector3<double>{0, 0, 1}, 0, 2
        ));
        return protein;
    }
}

TEST_CASE("GPU calculations agree with the CPU") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    bool device = on(true, [] () {return gpu::GPULoader::available();});
    WARN((device ? "comparing the CPU against the GPU device " + gpu::GPULoader::device_name() : "no GPU backend available, so this checks the CPU fallback"));

    Molecule protein("tests/files/2epe.pdb");
    protein.generate_new_hydration();

    SECTION("HistogramManagerMT") {
        check_manager<hist::HistogramManagerMT<false, false>>(protein);
        check_manager<hist::HistogramManagerMT<true, false>>(protein);
    }

    SECTION("HistogramManagerMTFFAvg") {
        check_manager<hist::HistogramManagerMTFFAvg<false, false>>(protein);
        check_manager<hist::HistogramManagerMTFFAvg<true, false>>(protein);
    }

    SECTION("HistogramManagerMTFFGrid") {
        check_manager<hist::HistogramManagerMTFFGrid<false>>(protein);
    }

    SECTION("SymmetryManagerMT") {
        auto symmetric = symmetric_2epe(protein.get_waters());
        check_manager<hist::SymmetryManagerMT<false, false>>(symmetric);
        check_manager<hist::SymmetryManagerMT<true, false>>(symmetric);
    }

    SECTION("incremental updates") {
        auto choice = GENERATE(
            settings::hist::HistogramManagerChoice::PartialHistogramManagerMT,
            settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT
        );
        auto waters = protein.get_waters();
        auto cpu = symmetric_2epe(waters);
        auto gpu = symmetric_2epe(waters);
        cpu.set_histogram_manager(choice);
        gpu.set_histogram_manager(choice);

        // each edit is applied to both molecules, and their histograms compared after it
        std::vector<std::pair<std::string, std::function<void(Molecule&)>>> edits = {
            {"initial", [] (Molecule&) {}},
            {"translate a body", [] (Molecule& m) {m.get_body(1).translate({2, 0, 0});}},
            {"change a symmetry", [] (Molecule& m) {
                static_cast<symmetry::CyclicSymmetry*>(m.get_body(0).symmetry().get(0))->_repeat_relation.translation = {0, 20, 0};
            }},
            {"change an atom", [] (Molecule& m) {
                m.get_body(2).get_atom(0).weight() = 2;
                m.get_body(2).get_signaller()->modified_internal();
            }},
            {"remove the hydration", [] (Molecule& m) {m.clear_hydration();}},
        };
        for (const auto& [name, edit] : edits) {
            INFO(name);
            edit(cpu);
            edit(gpu);
            auto h_cpu = on(false, [&] () {return cpu.get_histogram();});
            auto h_gpu = on(true,  [&] () {return gpu.get_histogram();});
            check(*h_cpu, *h_gpu);
        }
    }
}
