#pragma once

#include <md/utility/files/all.h>

#include <tuple>

namespace ausaxs::md {
    struct MDRunResult {
        MDRunResult(const io::Folder& folder) {
            gro = GROFile(folder + "prod.gro");
            edr = EDRFile(folder + "prod.edr");
            xtc = XTCFile(folder + "prod.xtc");
        }

        operator std::tuple<GROFile, EDRFile, XTCFile>() const {
            return std::make_tuple(gro, edr, xtc);
        }

        GROFile gro;
        EDRFile edr;
        XTCFile xtc;
    };

    struct SAXSRunResult {
        SAXSRunResult(const io::Folder& folder) {
            dat = DATFile(folder + "saxs.dat");
            xvg = XVGFile(folder + "waxs_final.xvg");
        }

        operator std::tuple<DATFile, XVGFile>() const {
            return std::make_tuple(dat, xvg);
        }

        DATFile dat;
        XVGFile xvg;
    };
}