#pragma once

#include <md/utility/files/MDPFile.h>

#include <vector>

namespace md {
    namespace MDPOptions {
        namespace detail {
            struct OptionVal {
                OptionVal(std::string name, std::string value) : name(name), value(value) {}
                OptionVal(std::string name, double value) : name(name), value(std::to_string(value)) {}
                OptionVal(std::string name, int value) : name(name), value(std::to_string(value)) {}
                OptionVal(std::string name, unsigned int value) : name(name), value(std::to_string(value)) {}
                std::string name, value;
            };

            struct Option {
                Option(std::string name) : name(name) {}
                Option(const char* name) : name(name) {}

                OptionVal operator=(std::string value) const {return OptionVal(name, value);}
                OptionVal operator=(double value) const {return OptionVal(name, value);}
                OptionVal operator=(int value) const {return OptionVal(name, value);}
                OptionVal operator=(unsigned int value) const {return OptionVal(name, value);}

                std::string name;
            };
        }

        const inline static detail::Option include     = "include";                   // include files
        const inline static detail::Option define      = "define";                    // define symbols

        const inline static detail::Option integrator  = "integrator";                // integrator
        const inline static detail::Option tinit       = "tinit";                     // initial time
        const inline static detail::Option dt          = "dt";                        // time step
        const inline static detail::Option nsteps      = "nsteps";                    // number of steps
        const inline static detail::Option nstcomm     = "nstcomm";                   // center of mass motion removal frequency
        const inline static detail::Option comm_mode   = "comm_mode";                 // center of mass motion removal mode

        const inline static detail::Option emtol       = "emtol";                     // energy tolerance
        const inline static detail::Option emstep      = "emstep";                    // energy step

        const inline static detail::Option rtpi        = "rtpi";                      // test particle insertion

        const inline static detail::Option nstlog      = "nstlog";                    // log output frequency
        const inline static detail::Option nstenergy   = "nstenergy";                 // energy output frequency
        const inline static detail::Option nstcalcenergy = "nstcalcenergy";           // energy calculation frequency
        const inline static detail::Option nstxout_compressed = "nstxout-compressed"; // compressed trajectory output frequency

        const inline static detail::Option cutoff_scheme = "cutoff-scheme";           // neighbor list cut-off scheme
        const inline static detail::Option nstlist     = "nstlist";                   // neighbor list frequency
        const inline static detail::Option pbc         = "pbc";                       // periodic boundary conditions
        const inline static detail::Option rlist       = "rlist";                     // short-range neighbor list cut-off

        const inline static detail::Option coulombtype = "coulombtype";               // electrostatics method
        const inline static detail::Option rcoulomb    = "rcoulomb";                  // short-range electrostatic cut-off
        const inline static detail::Option coulomb_modifier = "coulomb-modifier";     // electrostatics cut-off modifier
        const inline static detail::Option vdw_type    = "vdwtype";                   // van der Waals method
        const inline static detail::Option rvdw        = "rvdw";                      // short-range van der Waals cut-off
        const inline static detail::Option vdw_modifier= "vdw-modifier";              // van der Waals cut-off modifier
        const inline static detail::Option dispcorr    = "dispcorr";                  // long-range dispersion correction

        const inline static detail::Option tcoupl      = "tcoupl";                    // temperature coupling
        const inline static detail::Option tc_grps     = "tc_grps";                   // temperature coupling groups
        const inline static detail::Option tau_t       = "tau_t";                     // temperature coupling time
        const inline static detail::Option ref_t       = "ref_t";                     // reference temperature

        const inline static detail::Option pcoupl      = "pcoupl";                    // pressure coupling
        const inline static detail::Option pcoupltype  = "pcoupltype";                // pressure coupling type
        const inline static detail::Option tau_p       = "tau_p";                     // pressure coupling time
        const inline static detail::Option ref_p       = "ref_p";                     // reference pressure
        const inline static detail::Option compressibility = "compressibility";       // isothermal compressibility
        const inline static detail::Option refcoord_scaling = "refcoord_scaling";     // reference coordinate scaling

        const inline static detail::Option gen_vel     = "gen_vel";                   // generate velocities
        const inline static detail::Option gen_temp    = "gen_temp";                  // temperature for velocity generation

        const inline static detail::Option constraints = "constraints";                   // constraints
        const inline static detail::Option constraint_algorithm = "constraint_algorithm"; // constraints
        const inline static detail::Option lincs_order = "lincs-order";                   // LINCS order
        const inline static detail::Option lincs_iter = "lincs-iter";                     // LINCS iterations

        const inline static detail::Option scatt_coupl      = "scatt-coupl";
        const inline static detail::Option waxs_solute      = "waxs-solute";
        const inline static detail::Option waxs_solvent     = "waxs-solvent";
        const inline static detail::Option waxs_solvdens    = "waxs-solvdens";
        const inline static detail::Option waxs_rotfit      = "waxs-rotfit";
        const inline static detail::Option waxs_pbcatom     = "waxs-pbcatom";
        const inline static detail::Option waxs_fc          = "waxs-fc";
        const inline static detail::Option waxs_weights     = "waxs-weights";
        const inline static detail::Option waxs_nstcalc     = "waxs-nstcalc";
        const inline static detail::Option waxs_nstlog      = "waxs-nstlog";
        const inline static detail::Option waxs_nq          = "waxs-nq";
        const inline static detail::Option waxs_startq      = "waxs-startq";
        const inline static detail::Option waxs_endq        = "waxs-endq";
        const inline static detail::Option waxs_nsphere     = "waxs-nsphere";
        const inline static detail::Option waxs_solvdens_uncert = "waxs-solvdens-uncert";
        const inline static detail::Option waxs_correct_buffer = "waxs-correct-buffer";
    }

    class MDPCreator {
        public:
            MDPCreator& add(const MDPOptions::detail::OptionVal& option);
            std::string get(const MDPOptions::detail::Option& option) const;
            MDPFile write(const std::string& name) const;

        private:
            std::vector<MDPOptions::detail::OptionVal> options;
    };

    /**
     * @brief Template for energy minimization MDP files.
     * 
     * Fully works out of the box, but can be modified by adding options.
     */
    struct EMMDPCreator : public MDPCreator {EMMDPCreator();};

    /**
     * @brief Template for equilibration MDP files.
     */
    struct EQMDPCreator : public MDPCreator {EQMDPCreator();};

    /**
     * @brief Template for equilibration MDP files for molecules.
     * 
     * Fully works out of the box, but can be modified by adding options.
     */
    struct EQMDPCreatorMol : public EQMDPCreator {EQMDPCreatorMol();};

    /**
     * @brief Template for equilibration MDP files for solvents.
     * 
     * Fully works out of the box, but can be modified by adding options.
     */
    struct EQMDPCreatorSol : public EQMDPCreator {EQMDPCreatorSol();};

    /**
     * @brief Template for production MDP files.
     */
    struct PRMDPCreator : public MDPCreator {PRMDPCreator();};

    /**
     * @brief Template for production MDP files for molecules.
     * 
     * Fully works out of the box, but can be modified by adding options.
     */
    struct PRMDPCreatorMol : public PRMDPCreator {PRMDPCreatorMol();};

    /**
     * @brief Template for production MDP files for solvents.
     * 
     * Fully works out of the box, but can be modified by adding options.
     */
    struct PRMDPCreatorSol : public PRMDPCreator {PRMDPCreatorSol();};

    /**
     * @brief Template for SAXS MDP files.
     */
    struct SAXSMDPCreator : public MDPCreator {SAXSMDPCreator();};

    /**
     * @brief Template for SAXS MDP files for molecules.
     * 
     * The options
     *  - waxs-pbcatom
     *  - waxs-nsphere
     * must manually be set.
     */
    struct SAXSMDPCreatorMol : public SAXSMDPCreator {SAXSMDPCreatorMol();};

    /**
     * @brief Template for SAXS MDP files for solvents.
     * 
     * Fully works out of the box, but can be modified by adding options.
     */
    struct SAXSMDPCreatorSol : public SAXSMDPCreator {SAXSMDPCreatorSol();};

    struct TimeAnalysisMDPCreator : public MDPCreator {TimeAnalysisMDPCreator();};
    struct TimeAnalysisMDPCreatorMol : public TimeAnalysisMDPCreator {TimeAnalysisMDPCreatorMol();};
    struct TimeAnalysisMDPCreatorSol : public TimeAnalysisMDPCreator {TimeAnalysisMDPCreatorSol();};

    struct FrameAnalysisMDPCreator : public MDPCreator {FrameAnalysisMDPCreator();};
    struct FrameAnalysisMDPCreatorMol : public FrameAnalysisMDPCreator {FrameAnalysisMDPCreatorMol();};
    struct FrameAnalysisMDPCreatorSol : public FrameAnalysisMDPCreator {FrameAnalysisMDPCreatorSol();};
}