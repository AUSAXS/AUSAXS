#pragma once

#include <md/programs/options/water_models/IWaterModel.h>

namespace ausaxs::md::option::water_model {
    struct TIP4P : IWaterModel {
        std::string filename() const override { return "tip4p"; }
        std::string name() const override { return "TIP4P"; }
        std::string info() const override { return "TIP 4-point"; }
        std::string get_gro_file_content() const override {
            const char* gro_content = 
                #include "md/programs/options/water_models/tip4p.gro"
            ;
            return std::string(gro_content);
        }

        std::string get_itp_file_content() const override {
            return 
R"([ moleculetype ]
SOL		2

[ atoms ]
  1   OW_tip4p    1       SOL       OW       1       0        16.00000
  2   HW_tip4p    1       SOL       HW1      1       0.52      1.00800
  3   HW_tip4p    1       SOL       HW2      1       0.52      1.00800
  4   MW          1       SOL       MW       1      -1.04      0.00000

#ifndef FLEXIBLE

[ settles ]
1	1	0.09572	0.15139

#else

[ bonds ]
1       2       1       0.09572 502416.0 0.09572        502416.0 
1       3       1       0.09572 502416.0 0.09572        502416.0 
        
[ angles ]
2       1       3       1       104.52  628.02  104.52  628.02  

#endif

[ virtual_sites3 ]
4       1       2       3       1       0.128012065     0.128012065

[ exclusions ]
1	2	3	4
2	1	3	4
3	1	2	4
4	1	2	3

#ifdef SCATTER
[ scattering_params ]
1      1     CROMER_MANN_Owat
2      1     CROMER_MANN_Hwat
3      1     CROMER_MANN_Hwat

1      2     NEUTRON_SCATT_LEN_O
2      2     NSL_H_DEUTERATABLE
3      2     NSL_H_DEUTERATABLE
#endif)";
        }
    };
}