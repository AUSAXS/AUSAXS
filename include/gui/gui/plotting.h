#pragma once

#include <io/ExistingFile.h>

#include <array>
#include <fstream>

namespace resources {
    extern const std::array<unsigned char, 3759> plot_py;
    extern const std::array<unsigned char, 16726> plot_helper_py;

    inline io::ExistingFile generate_plotting_script() {
        io::File file("resources/plot.py");
        if (file.exists()) {return file;}
        else {file.create();}

        // plot.py
        std::ofstream out(file.path(), std::ios::binary);
        out.write(reinterpret_cast<const char*>(plot_py.data()), plot_py.size());

        // plot_helper.py
        std::ofstream helper("resources/plot_helper.py");
        helper.write(reinterpret_cast<const char*>(plot_helper_py.data()), plot_helper_py.size());
        return file;
    }
}