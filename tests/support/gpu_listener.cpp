#include <catch2/catch_test_case_info.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

#include <gpu/GPULoader.h>
#include <settings/GeneralSettings.h>

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <string_view>

using namespace ausaxs;

/**
 * Reruns the tests tagged [gpu] on the GPU backend: when the environment variable AUSAXS_TEST_GPU is set, the GPU is turned on before each
 * of them starts, so that every calculator they construct dispatches to the backend. Linked into every test executable.
 *
 * The GPU code falls back to the CPU whenever the backend cannot be used, and the tests would then pass without testing it. So a requested
 * GPU run fails outright if the backend is not usable when a test starts, or has been taken out of use by a failure when it ends.
 *
 * On a machine without a GPU, the backend can run on its host CPU device instead: `ausaxs gpu install --device host`.
 */
namespace {
    bool requested() {
        const char* value = std::getenv("AUSAXS_TEST_GPU");
        return value != nullptr && std::string_view(value) != "" && std::string_view(value) != "0";
    }

    bool tagged(const Catch::TestCaseInfo& info) {
        for (const auto& tag : info.tags) {
            if (tag.original == "gpu") {return true;}
        }
        return false;
    }

    class GPUListener : public Catch::EventListenerBase {
        public:
            using EventListenerBase::EventListenerBase;

            void testCaseStarting(const Catch::TestCaseInfo& info) override {
                if (!requested() || !tagged(info)) {return;}
                settings::general::gpu = true;
                if (!gpu::GPULoader::available()) {
                    throw std::runtime_error("AUSAXS_TEST_GPU is set, but the GPU backend cannot be used: " + gpu::GPULoader::device_name());
                }
            }

            void testCaseEnded(const Catch::TestCaseStats& stats) override {
                if (!requested() || !tagged(*stats.testInfo)) {return;}
                if (!gpu::GPULoader::available()) {
                    throw std::runtime_error("the GPU backend failed during \"" + stats.testInfo->name + "\", which therefore ran on the CPU");
                }
            }
    };
}

CATCH_REGISTER_LISTENER(GPUListener)
