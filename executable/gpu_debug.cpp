#include <webgpu/webgpu.hpp>
#include <gpu/WebGPU/WebGPUController.h>

int main(int, char const *[]) {
    std::vector<Atom> atoms = {
        Atom({-1, -1, -1}, 1), Atom({-1, 1, -1}, 1),
        Atom({ 1, -1, -1}, 1), Atom({ 1, 1, -1}, 1),
        Atom({-1, -1,  1}, 1), Atom({-1, 1,  1}, 1),
        Atom({ 1, -1,  1}, 1), Atom({ 1, 1,  1}, 1)
    };

    ausaxs::gpu::WebGPU<true> webgpu;
    webgpu.submit_self(atoms);
    auto res = webgpu.run();
    std::cout << "results: " << std::endl;
    std::cout << "\tself: " << res.self.size() << std::endl;
    std::cout << "\tcross: " << res.cross.size() << std::endl;

    auto res1 = res.self.at(0);
    std::cout << "Results:" << std::endl;
    for (int i = 0; i < 20; ++i) {
        std::cout << res1[i] << std::endl;
    }

    return 0;
}