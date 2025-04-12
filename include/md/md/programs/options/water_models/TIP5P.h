#pragma once

#include <md/programs/options/water_models/IWaterModel.h>

namespace ausaxs::md::option::water_model {
    struct TIP5P : IWaterModel {
        std::string filename() const override { return "tip5p"; }
        std::string name() const override { return "TIP5P"; }
        std::string info() const override { return "TIP 5-point (see https://gitlab.com/gromacs/gromacs/-/issues/1348 for issues)"; }
        std::string get_gro_file_content() const override {
            const char* gro_content = 
                #include "md/programs/options/water_models/tip5p.gro"
            ;
            return std::string(gro_content);
        }
        std::string get_itp_file_content() const override {
            return 
R"(
[ moleculetype ]
SOL		2

[ atoms ]
  1   OW_tip5p    1       SOL       OW       1       0        16.00000
  2   HW_tip5p    1       SOL       HW1      1       0.241     1.00800
  3   HW_tip5p    1       SOL       HW2      1       0.241     1.00800
  4   MW          1       SOL       LP1      1      -0.241     0.00000
  5   MW          1       SOL       LP2      1      -0.241     0.00000

#ifndef FLEXIBLE
[ settles ]
1	1	0.09572	0.15139

#else
[ bonds ]
; i     j       funct   length  force.c.
1       2       1       0.09572 502416.0 0.09572        502416.0 
1       3       1       0.09572 502416.0 0.09572        502416.0 
        
[ angles ]
2       1       3       1       104.52  628.02  104.52  628.02  
#endif

[ virtual_sites3 ]
4      1       2       3       4        -0.344908262    -0.34490826     -6.4437903493
5      1       2       3       4        -0.344908262    -0.34490826     6.4437903493

[ exclusions ]
1	2	3	4	5     
2	1	3	4	5
3	1	2	4	5
4	1	2	3	5
5	1	2	3	4

#ifdef SCATTER
[ scattering_params ]
1	1	CROMER_MANN_Owat
2	1	CROMER_MANN_Hwat
3	1	CROMER_MANN_Hwat
#endif
)";
        }
    };
}