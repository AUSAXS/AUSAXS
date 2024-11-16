#pragma once

#include <md/shell/Jobscript.h>
#include <md/utility/files/all.h>

#include <string>

namespace ausaxs::md {
    struct MDRunResult {
        MDRunResult(const std::string& folder) {
            if (folder.back() == '/') {
                gro = GROFile(folder + "prod.gro");
                edr = EDRFile(folder + "prod.edr");
                xtc = XTCFile(folder + "prod.xtc");
            } else {
                gro = GROFile(folder + "/prod.gro");
                edr = EDRFile(folder + "/prod.edr");
                xtc = XTCFile(folder + "/prod.xtc");
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