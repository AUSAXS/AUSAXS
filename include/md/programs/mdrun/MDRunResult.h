#pragma once

#include <shell/Jobscript.h>
#include <utility/files/all.h>

#include <string>
#include <memory>

namespace gmx {
    struct MDRunResult {
        MDRunResult(const std::string& folder) {
            if (folder.back() == '/') {
                gro = GROFile(folder + "confout.gro");
                edr = EDRFile(folder + "ener.edr");
                xtc = XTCFile(folder + "traj_comp.xtc");
            } else {
                gro = GROFile(folder + "/confout.gro");
                edr = EDRFile(folder + "/ener.edr");
                xtc = XTCFile(folder + "/traj_comp.xtc");
            }
        }
        GROFile gro;
        EDRFile edr;
        XTCFile xtc;

        operator std::tuple<GROFile, EDRFile, XTCFile>() const {
            return std::make_tuple(gro, edr, xtc);
        }
    };

    struct SAXSRunResult {
        SAXSRunResult(const std::string& folder) {
            if (folder.back() == '/') {
                dat = DATFile(folder + "saxs.dat");
                xvg = XVGFile(folder + "waxs_final.xvg");
            } else {
                dat = DATFile(folder + "/saxs.dat");
                xvg = XVGFile(folder + "/waxs_final.xvg");
            }
        }
        DATFile dat;
        XVGFile xvg;

        operator std::tuple<DATFile, XVGFile>() const {
            return std::make_tuple(dat, xvg);
        }
    };
}