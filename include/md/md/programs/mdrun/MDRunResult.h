#pragma once

#include <md/utility/files/all.h>
#include <io/Folder.h>

#include <tuple>
#include <vector>

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

    /**
     * @brief Result of a multi-replica mdrun (e.g. PLUMED -multidir).
     *
     * Constructed from the shared base folder that contains one replica
     * subdirectory per replica (named rep0, rep1, …). Each subdirectory is
     * expected to hold the standard prod.gro / prod.edr / prod.xtc output.
     */
    struct MultiMDRunResult {
        MultiMDRunResult(const io::Folder& base_folder) {
            for (auto& dir : base_folder.directories()) {
                std::string p = dir.path();
                auto slash = p.rfind('/');
                auto dirname = (slash == std::string::npos) ? p : p.substr(slash + 1);
                if (dirname.size() >= 3 && dirname.substr(0, 3) == "rep") {
                    replicas.emplace_back(dir);
                }
            }
        }

        std::vector<MDRunResult> replicas;
    };

    struct SAXSRunResult {
        SAXSRunResult(const io::Folder& folder) {
            // dat = DATFile(folder + "saxs.dat");
            xvg = XVGFile(folder + "prod_final.xvg");
        }

        operator std::tuple<DATFile, XVGFile>() const {
            return std::make_tuple(dat, xvg);
        }

        DATFile dat;
        XVGFile xvg;
    };
}