#include <md/programs/options/forcefields/IForcefield.h>
#include <md/programs/options/forcefields/AMBER99SB.h>
#include <md/programs/options/forcefields/AMBER99SB_ILDN.h>
#include <io/Folder.h>
#include <utility/Exceptions.h>
#include <settings/MDSettings.h>

using namespace ausaxs::md::option;

std::unique_ptr<IForcefield> IForcefield::construct(Forcefield ff) {
    switch (ff) {
        case Forcefield::AMBER99SB:
            return std::make_unique<force_field::AMBER99SB>();
        case Forcefield::AMBER99SB_ILDN:
            return std::make_unique<force_field::AMBER99SB_ILDN>();
        default:
            throw except::unknown_argument("gmx::to_string: Unknown force field. Did you forget to add it to the enum?");
    }
}

bool IForcefield::exists() const {
    io::Folder f(settings::md::gmx_top_path() + filename() + ".ff");
    return f.exists();
}
