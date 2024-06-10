#include <md/utility/files/MDPCreator.h>
#include <md/utility/Exceptions.h>

#include <fstream>
#include <iomanip>

using namespace md;

MDPCreator& MDPCreator::add(const MDPOptions::detail::OptionVal& option) {
    options.push_back(option);
    return *this;
}

MDPFile MDPCreator::write(const std::string& name) const {
    MDPFile file(name);
    std::ofstream out(file.path);                
    for (const auto& option : options) {
        out << std::left << std::setw(30) << option.name << " = " << option.value << std::endl;
    }
    return file;
}

std::string MDPCreator::get(const MDPOptions::detail::Option& option) const {
    for (const auto& o : options) {
        if (o.name == option.name) {
            return o.value;
        }
    }
    throw except::io_error("MDPCreator::get: option \"" + option.name + "\" not found.");
}

EMMDPCreator::EMMDPCreator() {
    add(MDPOptions::define = "-DFLEXIBLE");
    add(MDPOptions::integrator = "steep");
    add(MDPOptions::dt = 0.002);
    add(MDPOptions::nsteps = 100);
    add(MDPOptions::emtol = 1e-6);
    add(MDPOptions::emstep = 0.001);
    add(MDPOptions::lincs_order = 8);
    add(MDPOptions::lincs_iter = 2);
}

EQMDPCreator::EQMDPCreator() {
    add(MDPOptions::integrator = "md");
    add(MDPOptions::dt = 0.002);
    add(MDPOptions::nsteps = 50000);
    add(MDPOptions::nstxout_compressed = 5000);
    add(MDPOptions::cutoff_scheme = "Verlet");
    add(MDPOptions::nstlist = 20);
    add(MDPOptions::coulombtype = "PME");
    add(MDPOptions::coulomb_modifier = "Potential-shift");
    add(MDPOptions::vdw_type = "Cut-off");
    add(MDPOptions::vdw_modifier = "Potential-shift");
    add(MDPOptions::dispcorr = "EnerPres");
    add(MDPOptions::tcoupl = "v-rescale");
    add(MDPOptions::pcoupl = "Berendsen");
    add(MDPOptions::pcoupltype = "isotropic");
    add(MDPOptions::tau_p = 1.0);
    add(MDPOptions::ref_p = 1.0);
    add(MDPOptions::compressibility = 4.5e-5);
    add(MDPOptions::refcoord_scaling = "com");
    add(MDPOptions::gen_vel = "yes");
    add(MDPOptions::constraints = "h-bonds");
}

EQMDPCreatorMol::EQMDPCreatorMol() : EQMDPCreator() {
    add(MDPOptions::define = "-DPOSRES");
    add(MDPOptions::tc_grps = "Protein Water_and_ions");
    add(MDPOptions::tau_t = "0.1 0.1");
    add(MDPOptions::ref_t = "298.15 298.15");
}

EQMDPCreatorSol::EQMDPCreatorSol() : EQMDPCreator() {
    add(MDPOptions::tc_grps = "Water");
    add(MDPOptions::tau_t = "0.1");
    add(MDPOptions::ref_t = "298.15");
}

PRMDPCreator::PRMDPCreator() {
    add(MDPOptions::integrator = "md");
    add(MDPOptions::dt = 0.002);
    add(MDPOptions::nsteps = 500000);
    add(MDPOptions::nstxout_compressed = 1000);
    add(MDPOptions::cutoff_scheme = "Verlet");
    add(MDPOptions::nstlist = 20);
    add(MDPOptions::coulombtype = "PME");
    add(MDPOptions::coulomb_modifier = "Potential-shift-Verlet");
    add(MDPOptions::vdw_type = "Cut-off");
    add(MDPOptions::vdw_modifier = "Potential-shift-Verlet");
    add(MDPOptions::dispcorr = "EnerPres");
    add(MDPOptions::tcoupl = "v-rescale");
    add(MDPOptions::pcoupl = "Parrinello-Rahman");
    add(MDPOptions::pcoupltype = "isotropic");
    add(MDPOptions::tau_p = 5.0);
    add(MDPOptions::ref_p = 1.0);
    add(MDPOptions::compressibility = 4.5e-5);
    add(MDPOptions::refcoord_scaling = "com");
    add(MDPOptions::gen_vel = "yes");
    add(MDPOptions::constraints = "h-bonds");
}

PRMDPCreatorMol::PRMDPCreatorMol() : PRMDPCreator() {
    add(MDPOptions::define = "-DPOSRESBACKBONE");
    add(MDPOptions::tc_grps = "Protein Water_and_ions");
    add(MDPOptions::tau_t = "0.1 0.1");
    add(MDPOptions::ref_t = "298.15 298.15");
}

PRMDPCreatorSol::PRMDPCreatorSol() : PRMDPCreator() {
    add(MDPOptions::tc_grps = "Water");
    add(MDPOptions::tau_t = "0.1");
    add(MDPOptions::ref_t = "298.15");
}

SAXSMDPCreator::SAXSMDPCreator() {
    add(MDPOptions::define = "-DSCATTER");
    // add(MDPOptions::integrator = "md");
    // add(MDPOptions::dt = 0.002);
    // add(MDPOptions::nsteps = 500000);
    // add(MDPOptions::nstlog = 5000);
    // add(MDPOptions::nstenergy = 5000);
    // add(MDPOptions::nstcomm = 1);
    // add(MDPOptions::nstxout_compressed = 5000);

    // add(MDPOptions::cutoff_scheme = "Verlet");
    // add(MDPOptions::nstlist = 20);

    // add(MDPOptions::coulombtype = "PME");
    // add(MDPOptions::coulomb_modifier = "Potential-shift-Verlet");
    // add(MDPOptions::vdw_type = "Cut-off");
    // add(MDPOptions::vdw_modifier = "Potential-shift-Verlet");
    // add(MDPOptions::dispcorr = "EnerPres");

    // add(MDPOptions::tcoupl = "v-rescale");
    // add(MDPOptions::pcoupl = "Berendsen");
    // add(MDPOptions::pcoupltype = "isotropic");
    // add(MDPOptions::ref_p = 1.0);
    // add(MDPOptions::tau_p = 1.0);
    // add(MDPOptions::compressibility = 4.5e-5);
    // add(MDPOptions::refcoord_scaling = "com");

    // add(MDPOptions::constraints = "h-bonds");
    // add(MDPOptions::constraint_algorithm = "Lincs");

    add(MDPOptions::scatt_coupl = "xray");
    add(MDPOptions::waxs_solvdens = 334);
    add(MDPOptions::waxs_fc = 0);
    add(MDPOptions::waxs_weights = "exp+calc+solvdens");
    add(MDPOptions::waxs_nstcalc = 125);
    add(MDPOptions::waxs_nstlog = 5000);
    add(MDPOptions::waxs_nq = 101);
    add(MDPOptions::waxs_startq = 0);
    add(MDPOptions::waxs_endq = 5);
    add(MDPOptions::waxs_solvdens_uncert = 0.005);
}

SAXSMDPCreatorMol::SAXSMDPCreatorMol() : SAXSMDPCreator() {
    add(MDPOptions::tc_grps = "Protein Water_and_Ions");
    add(MDPOptions::ref_t = "298.15 298.15");
    add(MDPOptions::tau_t = "0.1 0.1");
    add(MDPOptions::waxs_solute = "Prot-Masses");
    add(MDPOptions::waxs_solvent = "RealWater_and_Ions");
    add(MDPOptions::waxs_rotfit = "Prot-Masses");
}

SAXSMDPCreatorSol::SAXSMDPCreatorSol() : SAXSMDPCreator() {
    add(MDPOptions::tc_grps = "Water");
    add(MDPOptions::tau_t = "0.1");
    add(MDPOptions::ref_t = "298.15");
    add(MDPOptions::waxs_solute = "");
    add(MDPOptions::waxs_solvent = "RealWater");
    add(MDPOptions::waxs_rotfit = "");
    add(MDPOptions::waxs_pbcatom = 0);
    add(MDPOptions::waxs_nsphere = 1200);
    add(MDPOptions::waxs_correct_buffer = "yes");
}

TimeAnalysisMDPCreator::TimeAnalysisMDPCreator() {
    add(MDPOptions::integrator = "md");
    add(MDPOptions::dt = 0.002);
    add(MDPOptions::nsteps = 5000000);
    add(MDPOptions::nstxout_compressed = 1000);
    add(MDPOptions::cutoff_scheme = "Verlet");
    add(MDPOptions::nstlist = 20);
    add(MDPOptions::coulombtype = "PME");
    add(MDPOptions::coulomb_modifier = "Potential-shift-Verlet");
    add(MDPOptions::vdw_type = "Cut-off");
    add(MDPOptions::vdw_modifier = "Potential-shift-Verlet");
    add(MDPOptions::dispcorr = "EnerPres");
    add(MDPOptions::tcoupl = "v-rescale");
    add(MDPOptions::pcoupl = "Parrinello-Rahman");
    add(MDPOptions::pcoupltype = "isotropic");
    add(MDPOptions::tau_p = 5.0);
    add(MDPOptions::ref_p = 1.0);
    add(MDPOptions::compressibility = 4.5e-5);
    add(MDPOptions::refcoord_scaling = "com");
    add(MDPOptions::gen_vel = "yes");
    add(MDPOptions::constraints = "h-bonds");
}

TimeAnalysisMDPCreatorMol::TimeAnalysisMDPCreatorMol() {
    add(MDPOptions::define = "-DPOSRESBACKBONE");
    add(MDPOptions::tc_grps = "Protein Water_and_Ions");
    add(MDPOptions::tau_t = "0.1 0.1");
    add(MDPOptions::ref_t = "298.15 298.15");
}

TimeAnalysisMDPCreatorSol::TimeAnalysisMDPCreatorSol() {
    add(MDPOptions::tc_grps = "Water");
    add(MDPOptions::tau_t = "0.1");
    add(MDPOptions::ref_t = "298.15");
}

FrameAnalysisMDPCreator::FrameAnalysisMDPCreator() {
    add(MDPOptions::integrator = "md");
    add(MDPOptions::dt = 0.002);
    add(MDPOptions::nsteps = 500000);
    add(MDPOptions::nstxout_compressed = 250);
    add(MDPOptions::cutoff_scheme = "Verlet");
    add(MDPOptions::nstlist = 20);
    add(MDPOptions::coulombtype = "PME");
    add(MDPOptions::coulomb_modifier = "Potential-shift-Verlet");
    add(MDPOptions::vdw_type = "Cut-off");
    add(MDPOptions::vdw_modifier = "Potential-shift-Verlet");
    add(MDPOptions::dispcorr = "EnerPres");
    add(MDPOptions::tcoupl = "v-rescale");
    add(MDPOptions::pcoupl = "Parrinello-Rahman");
    add(MDPOptions::pcoupltype = "isotropic");
    add(MDPOptions::tau_p = 5.0);
    add(MDPOptions::ref_p = 1.0);
    add(MDPOptions::compressibility = 4.5e-5);
    add(MDPOptions::refcoord_scaling = "com");
    add(MDPOptions::gen_vel = "yes");
    add(MDPOptions::constraints = "h-bonds");    
}

FrameAnalysisMDPCreatorMol::FrameAnalysisMDPCreatorMol() {
    add(MDPOptions::define = "-DPOSRESBACKBONE");
    add(MDPOptions::tc_grps = "Protein Water_and_ions");
    add(MDPOptions::tau_t = "0.1 0.1");
    add(MDPOptions::ref_t = "298.15 298.15");
}

FrameAnalysisMDPCreatorSol::FrameAnalysisMDPCreatorSol() {
    add(MDPOptions::tc_grps = "Water");
    add(MDPOptions::tau_t = "0.1");
    add(MDPOptions::ref_t = "298.15");
}

