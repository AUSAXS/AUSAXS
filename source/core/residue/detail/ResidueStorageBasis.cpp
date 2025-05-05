#include <settings/GeneralSettings.h>

#include <string>
#include <fstream>

namespace ausaxs::residue::detail {
    void write_master_basis();
}

void ausaxs::residue::detail::write_master_basis() {
    std::string basis = R"(
#
LYS
O oc2 0
O oxt 1
N nz 3
C ce 2
O oc1 1
C cd 2
C cg 2
O o 0
C c 0
C ca 1
C cb 2
N n 1

#
VAL
O oc2 0
O oc1 1
O oxt 1
C cg2 3
C cg1 3
C cg 3
O o 0
C c 0
C ca 1
C cb 1
N n 1

#
PHE
O oc2 0
O oxt 1
N n 1
C cb 2
C ce1 1
C ca 1
O oc1 1
C c 0
O o 0
C cg 0
C cd 1
C ce2 1
C cd1 1
C cd2 1
C ce 1
C cz 1

#
GLY
O oc2 0
O oc1 1
O oxt 1
O o 0
C c 0
C ca 2
N n 1

#
ARG
O oc2 0
N n 1
C cb 2
C ca 1
O o 0
C cg 2
C cd 2
C cz 0
O oc1 1
N nh1 2
C c 0
N nh 2
N nh2 2
N ne 1
O oxt 1

#
CYS
O oc2 0
O oc1 1
O oxt 1
S sg 1
O o 0
C c 0
C ca 1
C cb 2
N n 1

#
GLU
O oc2 0
O oxt 1
O oe2 1
O oe1 0
O oe 0
O oc1 1
C cd 0
C cg 2
O o 0
C c 0
C ca 1
C cb 2
N n 1

#
LEU
O oc2 0
O oxt 1
C cd2 3
C cd1 3
O oc1 1
C cd 3
C cg 1
O o 0
C c 0
C ca 1
C cb 2
N n 1

#
ALA
O oc2 0
O oc1 1
O oxt 1
O o 0
C c 0
C ca 1
C cb 3
N n 1

#
HIS
O oc2 0
N n 1
C cb 2
C ce1 1
C ca 1
O o 0
C cg 0
N nd 1
O oc1 1
C c 0
N ne2 1
N nd1 1
C cd2 1
C ce 1
O oxt 1

#
ASP
O oc2 0
O oxt 1
O oc1 1
O od2 1
O od1 0
O od 0
C cg 0
O o 0
C c 0
C ca 1
C cb 2
N n 1

#
ASN
O oc2 0
O oxt 1
N nd2 2
O oc1 1
O od1 0
O od 0
C cg 0
O o 0
C c 0
C ca 1
C cb 2
N n 1

#
TYR
O oc2 0
O oxt 1
N n 1
C cb 2
C ce1 1
O oh 1
C ca 1
O oc1 1
C c 0
O o 0
C cg 0
C cd 1
C ce2 1
C cd1 1
C cd2 1
C ce 1
C cz 0

#
SER
O oc2 0
O oc1 1
O oxt 1
O og 1
O o 0
C c 0
C ca 1
C cb 2
N n 1

#
TRP
O oc2 0
C cz2 1
N n 1
C cb 2
C ca 1
N ne1 1
O oc1 1
C c 0
O o 0
C ch2 1
C ce3 1
C cg 0
C cd 1
C cz3 1
C ce2 0
C cd1 1
C cd2 0
O oxt 1
N ne 1

#
THR
O oc2 0
O oc1 1
O oxt 1
C cg2 3
O og 1
O o 0
O og1 1
C c 0
C ca 1
C cb 1
N n 1

#
GLN
O oc2 0
O oxt 1
O oe1 0
O oe 0
O oc1 1
N ne2 2
C cd 0
C cg 2
O o 0
C c 0
C ca 1
C cb 2
N n 1

#
ILE
O oc2 0
O oxt 1
C cd1 3
O oc1 1
C cd 3
C cg2 3
C cg1 2
C cg 2
O o 0
C c 0
C ca 1
C cb 1
N n 1

#
PRO
O oc2 0
O oxt 1
O oc1 1
C cd 2
C cg 2
O o 0
C c 0
C ca 1
C cb 2
N n 0

#
MET
O oc1 1
O oxt 1
C ce 3
O oc2 0
S sd 0
C cg 2
O o 0
C c 0
C ca 1
C cb 2
N n 1

# fake residue to denote Cl ions, used by YASARA (and by extension, WAXSiS)
CIM
cl cl 0

# same as above, but for Na ions
CIP
na na 0
)";

    std::string path = settings::general::residue_folder;
    std::ofstream file(path + "master.dat");
    file << basis;
    file.close();
}