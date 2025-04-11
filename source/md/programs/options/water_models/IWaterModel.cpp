#include <md/programs/options/water_models/IWaterModel.h>
#include <md/programs/options/water_models/TIP4P.h>
#include <md/programs/options/water_models/TIP4P2005.h>
#include <md/programs/options/water_models/TIP5P.h>
#include <io/File.h>
#include <utility/Exceptions.h>
#include <settings/MDSettings.h>

#include <fstream>
#include <cassert>

using namespace ausaxs;
using namespace ausaxs::md;

std::unique_ptr<option::IWaterModel> option::IWaterModel::construct(WaterModel wm) {
    switch (wm) {
        case WaterModel::TIP4P:
            return std::make_unique<option::water_model::TIP4P>();
        case WaterModel::TIP4P2005:
            return std::make_unique<water_model::TIP4P2005>();
        case WaterModel::TIP5P:
            return std::make_unique<water_model::TIP5P>();
        default:
            throw except::unknown_argument("gmx::to_string: Unknown water model. Did you forget to add it to the enum?");
    }
}

void option::IWaterModel::create(observer_ptr<option::IForcefield> ff) const {
    // create itp file
    io::File f(settings::md::gmx_top_path() + ff->filename() + ".ff/" + filename() + ".itp");
    assert(!f.exists() && "gmx::create: File already exists: " + f.path());
    f.create(get_file_content());

    // update watermodels.dat
    io::File dat(settings::md::gmx_top_path() + ff->filename() + ".ff/watermodels.dat");
    if (!dat.exists()) {
        throw except::io_error("gmx::create: File does not exist: " + dat.path());
    }

    // check if present in the file
    std::ifstream dat_file(dat.path());
    std::string line;
    while (std::getline(dat_file, line)) {
        if (line.starts_with(filename())) {
            return;
        }
    }
    dat_file.close();

    // not present in the file, append it
    std::ofstream dat_file_out(dat.path(), std::ios_base::app);
    dat_file_out << filename() << " " << name() << " " << info() << std::endl;
    dat_file_out.close();
    if (!dat_file_out.good()) {
        throw except::io_error("gmx::create: Could not write to file: " + dat.path());
    }
}

bool option::IWaterModel::exists(observer_ptr<option::IForcefield> ff) const {
    io::File f(settings::md::gmx_top_path() + ff->filename() + ".ff/" + filename() + ".itp");
    return f.exists();
}