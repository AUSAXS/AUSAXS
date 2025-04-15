#include "utility/Console.h"
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

void option::IWaterModel::create_gro() const {
    console::print_text("Creating gro file for water model " + filename() + "...");
    io::File gro(settings::md::gmx_top_path() + filename() + ".gro");
    assert(!gro.exists() && "gmx::create: File already exists!");
    gro.create(get_gro_file_content());
}

void option::IWaterModel::create_itp(observer_ptr<option::IForcefield> ff) const {
    console::print_text("Creating itp file for water model " + filename() + "...");

    // create itp file
    io::File f(settings::md::gmx_top_path() + ff->filename() + ".ff/" + filename() + ".itp");
    assert(!f.exists() && "gmx::create: File already exists!");
    f.create(get_itp_file_content());

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

bool option::IWaterModel::itp_exists(observer_ptr<option::IForcefield> ff) const {
    io::File f(settings::md::gmx_top_path() + ff->filename() + ".ff/" + filename() + ".itp");
    return f.exists();
}

bool option::IWaterModel::gro_exists() const {
    io::File f(settings::md::gmx_top_path() + filename() + ".gro");
    return f.exists();
}

bool option::IWaterModel::atomtypes_defined(observer_ptr<option::IForcefield> ff) const {
    io::File f(settings::md::gmx_top_path() + ff->filename() + ".ff/ffnonbonded.itp");
    assert(f.exists() && "gmx::atomtypes_defined: File does not exist!");

    std::ifstream in(f.path());
    std::string line;
    while (std::getline(in, line)) {
        if (line.find(filename()) != std::string::npos) {
            return true;
        }
    }
    in.close();
    return false;
}

void option::IWaterModel::define_atomtypes(observer_ptr<option::IForcefield> ff) const {
    console::print_text("Defining atomtypes for water model " + filename() + "...");

    io::File f(settings::md::gmx_top_path() + ff->filename() + ".ff/ffnonbonded.itp");
    assert(f.exists() && "gmx::define_atomtypes: File does not exist!");

    std::ofstream out(f.path(), std::ios_base::app);
    if (!out.good()) {
        throw except::io_error("gmx::define_atomtypes: Could not write to file: " + f.path());
    }
    out << get_defined_atomtypes();
    out.close();
}

void option::IWaterModel::ensure_exists(observer_ptr<option::IForcefield> ff) {
    if (!itp_exists(ff)) {
        create_itp(ff);
    }
    if (!gro_exists()) {
        create_gro();
    }
    if (!atomtypes_defined(ff)) {
        define_atomtypes(ff);
    }
}