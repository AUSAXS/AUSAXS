// SPDX-License-Identifier: LGPL-3.0-or-later

#include <container/ThreadLocalWrapper.h>
#include <gpu/GPULoader.h>
#include <grid/detail/GridExcludedVolume.h>
#include <hist/detail/GridExvFFT.h>
#include <io/File.h>
#include <math/Vector3.h>
#include <settings/GeneralSettings.h>
#include <shell/Command.h>
#include <utility/Curl.h>
#include <utility/MultiThreading.h>

#include <stdexcept>
#include <string>

namespace {
    void require(bool condition, const char* message) {
        if (!condition) {throw std::runtime_error(message);}
    }
}

int main() {
    using namespace ausaxs;
    auto* pool = utility::multi_threading::get_global_pool();
    require(pool->get_thread_ids().empty(), "serial pool created worker threads");
    container::ThreadLocalWrapper<int> local(0);
    require(local.size() == 1, "serial thread-local storage must have one entry");
    pool->detach_task([&] {pool->detach_task([&] {++local.get();});});
    pool->detach_loop(0, 4, [&](int i) {local.get() += i;});
    int block_sum = 0;
    pool->detach_blocks(3, 5, [&](int first, int last) {block_sum = first + last;});
    pool->detach_blocks(4, 4, [&](int, int) {throw std::runtime_error("empty block executed");});
    pool->submit_blocks(0, 3, [&](int first, int last) {block_sum += last - first;}).wait();
    pool->wait();
    require(local.merge() == 7 && block_sum == 11, "serial task result is incorrect");
    try {
        pool->detach_task([] {throw std::runtime_error("task failure");});
        throw std::runtime_error("task exception was lost");
    } catch (const std::runtime_error& e) {
        require(std::string(e.what()) == "task failure", "incorrect task exception");
    }

    settings::general::threads = 8;
    require(settings::general::threads == 1, "serial build accepted multiple threads");
    settings::general::gpu = true;
    require(settings::general::gpu == false && !gpu::GPULoader::available(), "serial build accepted GPU use");
    require(!gpu::GPULoader::probe("backend").error.empty(), "GPU probe must report its failure");
    require(!curl::download("https://example.invalid", io::File("/unavailable")), "offline download succeeded");
    try {
        shell::Command("true").execute();
        throw std::runtime_error("shell command unexpectedly succeeded");
    } catch (const std::exception& e) {
        require(std::string(e.what()).find("unavailable") != std::string::npos, "shell failure was not reported");
    }

#ifdef POCKETFFT_AVAILABLE
    grid::exv::GridExcludedVolume exv;
    exv.spacing = 2;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                exv.interior.emplace_back(2*i, 2*j, 2*k);
                exv.interior_sites.emplace_back(i, j, k);
            }
        }
    }
    const auto lattice = hist::detail::lattice::self_correlation(exv, 1, 100);
    double count = 0;
    for (int i = 0; i < lattice.size(); ++i) {count += lattice.index(i).count;}
    require(count == 64*63, "incorrect FFT pair count");
    require(lattice.index(2).count == 2*3*(3*4*4), "incorrect FFT nearest-neighbor count");
#endif
}
