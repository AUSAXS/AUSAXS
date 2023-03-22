#pragma once

#include <vector>
#include <string>
#include <memory>
#include <thread>

#include <math/Vector3.h>
#include <utility/Limit3D.h>
#include <utility/Type.h>
#include <utility/Exceptions.h>

// A small container of the various settings. These should be set *before* their respective classes are instantiated. 
namespace setting {
    struct general {
        static constexpr char residue_folder[] = "temp/residues/";                  // Download location for all ligand files. Must be constexpr. 
        inline static bool verbose = true;                                          // Whether to print out extra information.
        inline static unsigned int threads = std::thread::hardware_concurrency();   // The number of threads to use for parallelization.
        inline static std::string output = "output/";                               // The output directory.
        inline static bool keep_hydrogens = false;                                  // Whether to keep bound hydrogens when reading a structure.

        struct detail {
            inline static unsigned int job_size = 200; // The number of atoms to process in each job.
        };
    };

    struct grid {
        enum class PlacementStrategy {AxesStrategy, RadialStrategy, JanStrategy};
        enum class CullingStrategy {CounterStrategy, OutlierStrategy, RandomStrategy};

        inline static PlacementStrategy placement_strategy = PlacementStrategy::RadialStrategy; // The choice of placement algorithm.
        inline static CullingStrategy culling_strategy = CullingStrategy::CounterStrategy;      // The choice of culling algorithm. 

        inline static double percent_water = 0.1; // The number of generated water molecules as a percent of the number of atoms. 
        inline static double ra = 2.4;            // Radius of protein atoms. 
        inline static double rh = 1.5;            // Radius of water molecules.
        inline static double ra_effective = 2.4;  // Effective radius of protein atoms. This is based on the volume the average atom effectively occupies. 
        inline static double width = 1;           // Width of each bin of the grid used to represent this protein.
        inline static double scaling = 0.25;      // The percent increase in grid size in all dimensions when the grid size is automatically deduced based on an input vector of atoms. 
        inline static bool cubic = false;         // Whether to generate a cubic grid. This is primarily intended for rigid body optimization, to ensure there's enough space for all possible conformations. 

        inline static Limit3D axes = Limit3D(-250, 250, -250, 250, -250, 250); // Default axes for the grid 

        struct placement {
            inline static double min_score = 0.1; // (0.5 + min_score) is the minimum percentage of radial lines which must not intersect anything to place a water molecule
        };
    };

    struct protein {
        inline static bool center = true;               // Decides if the structure will be centered at origo. 
        inline static bool use_effective_charge = true; // Decides whether the charge of the displaced water will be included. 
    };

    struct axes {
        inline static unsigned int max_distance = 2000; // The maximum distance in the p(r) function in Ångström. Should always be much larger than the actual maximum distance. 
        inline static double distance_bin_width = 1;    // The width of each bin for the scattering plots.
        inline static unsigned int bins = 1000;         // The number of bins for the scattering plots.
        inline static double qmin = 1e-4;               // Lower limit on the used q-values
        inline static double qmax = 0.5;                // Upper limit on the used q-values
        inline static unsigned int skip = 0;            // The number of points to skip from the top of the scattering curve.
    };

    struct fit {
        inline static unsigned int N = 100; // Number of points sampled when discretizing a model scattering curve
        inline static bool verbose = false; // Whether to print the fit progress to the console
    };

    struct rigidbody {
        enum class TransformationStrategyChoice {RigidTransform, SingleTransform};
        enum class ParameterGenerationStrategyChoice {Simple};
        enum class BodySelectStrategyChoice {RandomSelect, RandomConstraintSelect, SequentialSelect};

        inline static TransformationStrategyChoice tsc = TransformationStrategyChoice::RigidTransform;
        inline static ParameterGenerationStrategyChoice pgsc = ParameterGenerationStrategyChoice::Simple;
        inline static BodySelectStrategyChoice bssc = BodySelectStrategyChoice::RandomSelect;

        inline static unsigned int iterations = 1000;   // The number of iterations to run the rigid body optimization for.
        inline static double bond_distance = 3;         // The maximum distance in Ångström between two atoms that allows for a constraint. 

        struct detail {
            inline static std::vector<int> constraints; // The residue ids to place a constraint at.
        };
    };

    struct crystal {
        enum class MillerGenerationChoice {All, Fibonacci, Reduced};
        inline static MillerGenerationChoice mgc = MillerGenerationChoice::All; // The choice of Miller index generation algorithm.
        inline static unsigned int h = 100;                                     // The maximum Miller index along the x direction.
        inline static unsigned int k = 100;                                     // The maximum Miller index along the y direction.
        inline static unsigned int l = 100;                                     // The maximum Miller index along the z direction.
        struct reduced {
            inline static double hkl_limit = 4; // The maximum Miller length to use in the Reduced algorithm. h^2 + k^2 + l^2 <= hkl_limit^2
        };
    };

    struct em {
        enum class CullingStrategyChoice {NoStrategy, CounterStrategy};
        inline static CullingStrategyChoice csc = CullingStrategyChoice::CounterStrategy; // The choice of culling algorithm. 

        inline static unsigned int sample_frequency = 1; // How often a bin is sampled in any direction. 
        inline static double concentration = 2;          // The concentration in mg/mL used when calculating the absolute intensity scale for simulations.
        inline static unsigned int charge_levels = 100;  // The number of partial histograms to utilize.

        inline static bool hydrate = true;               // Whether to hydrate the protein in the EM algorithm.
        inline static unsigned int evals = 100;          // Base number of evaluations used in the EM fitter. 

        inline static bool save_pdb = true;              // Whether to save the final atomic structure as a PDB file.
        inline static Limit alpha_levels = {0.5, 8};     // The range of alpha-levels to search.

        inline static bool fixed_weights = false;        // Whether to use fixed or dynamic weights for the EM algorithm. Fixed weights means that all atoms will have the same weight of 1. 

        struct simulation {
            inline static bool noise = true; // Whether to generate noise for the simulations. 
        };
    };

    struct plot {
        inline static std::string format = "png";       // The output format.

        struct image {
            inline static std::vector<double> contour = {}; // The contour levels for the image plots.
        };

        struct em {
            // Whether to plot the evaluated chi2 points. 
            // Produces 2 plots; one of the full landscape and another of the area near the minimum. 
            // The number of points is roughly determined by setting::em::evals
            inline static bool additional_plots = true; 
            inline static bool landscape = true; 
        };
    };

    /**
     * @brief Read the settings from a file.
     */
    void read(std::string path);

    /**
     * @brief Write the settings to a file.
     */
    void write(std::string path);

    /**
     * @brief Check if a settings file exists in the given directory, and read it if so.
     * @return True if a settings file was found and read.
     */
    bool discover(std::string path);

    namespace detail {
        // Super class for SmartOptions to allow polymorphic vectors
        struct ISmartOption {
            ISmartOption(std::vector<std::string> aliases) : aliases(aliases) {}
            std::vector<std::string> aliases; // A list of strings which will be recognized as this option. The first string will be used for output.

            /**
             * @brief Assign a new value to the setting represented by this SmartOption.
             */
            virtual void set(std::vector<std::string>) const = 0;

            /**
             * @brief Get a string representation of the setting represented by this SmartOption.
             */
            virtual std::string get() const = 0;
        };

        // Option specifier. This class maps an option enum to its recognized aliases and value storage. 
        template<typename T>
        struct SmartOption : public ISmartOption {
            SmartOption(std::vector<std::string> aliases, T& setting) : ISmartOption(aliases), setting(setting) {}
            T& setting; // A pointer to the actual setting value. This is so we can access it later. 

            /**
             * @brief Assign a new value to the setting represented by this SmartOption. 
             *        The string will be cast into a proper type. This method must be specialized for all possible setting types.
             */
            void set(std::vector<std::string>) const override {throw except::unexpected("Settings::SmartOption::set: Not implemented for type \"" + type(setting) + "\".");}

            /**
             * @brief Get a string representation of the setting represented by this SmartOption.
             *        This method must be specialized for all possible setting types.
             */
            std::string get() const override {throw except::unexpected("Settings::SmartOption::set: Not implemented for type \"" + type(setting) + "\".");}
        };

        /**
         * @brief Helper function for creating shared pointers to SmartOptions. This is to automatically deduce template types,
         *        which is not possible with the default shared_ptr constructor.
         */
        template<typename T>
        std::shared_ptr<SmartOption<T>> make_shared(std::vector<std::string> aliases, T& setting) {
            return std::make_shared<SmartOption<T>>(aliases, setting);
        }

        // The following vector lists the settings which will be made available for both reading and writing in the settings file. 
        inline static const std::vector<std::shared_ptr<ISmartOption>> options = {
            // figures
            make_shared({"format"}, setting::plot::format),

            // grid
            make_shared({"percent-water"}, setting::grid::percent_water),
            make_shared({"ra"}, setting::grid::ra),
            make_shared({"rh"}, setting::grid::rh),
            make_shared({"ra-effective"}, setting::grid::ra_effective),
            make_shared({"width"}, setting::grid::width),
            make_shared({"scaling"}, setting::grid::scaling),
            make_shared({"cubic"}, setting::grid::cubic),

            // protein
            make_shared({"center"}, setting::protein::center),
            make_shared({"effectivecharge", "effective-charge", "use-effective-charge"}, setting::protein::use_effective_charge),

            // axes
            make_shared({"scattering-intensity-plot-binned-width"}, setting::axes::distance_bin_width),
            make_shared({"qlow", "qmin"}, setting::axes::qmin),
            make_shared({"qhigh", "qmax"}, setting::axes::qmax),
            make_shared({"skip", "skip-rows"}, setting::axes::skip),

            // fit
            make_shared({"N"}, setting::fit::N),

            // rigidbody
            make_shared({"rigidbody-iterations"}, setting::rigidbody::iterations),
            make_shared({"bond-distance"}, setting::rigidbody::bond_distance),
            make_shared({"constraints"}, setting::rigidbody::detail::constraints),

            // em
            make_shared({"sample-frequency"}, setting::em::sample_frequency),
            make_shared({"concentration"}, setting::em::concentration),
            make_shared({"charge-levels"}, setting::em::charge_levels),
            make_shared({"noise"}, setting::em::simulation::noise),
            make_shared({"contour"}, setting::plot::image::contour),
            make_shared({"fixed-weights"}, setting::em::fixed_weights),

            // plot

        };

        // declare template specializations
        template<> std::string SmartOption<std::vector<int>>::get() const;
        template<> std::string SmartOption<std::vector<std::string>>::get() const;
        template<> std::string SmartOption<std::vector<double>>::get() const;
        template<> std::string SmartOption<std::string>::get() const;
        template<> std::string SmartOption<double>::get() const;
        template<> std::string SmartOption<int>::get() const;
        template<> std::string SmartOption<unsigned int>::get() const;
        template<> std::string SmartOption<bool>::get() const;
        template<> void SmartOption<std::string>::set(std::vector<std::string> str) const;
        template<> void SmartOption<bool>::set(std::vector<std::string> str) const;
        template<> void SmartOption<double>::set(std::vector<std::string> str) const;
        template<> void SmartOption<int>::set(std::vector<std::string> str) const;
        template<> void SmartOption<unsigned int>::set(std::vector<std::string> str) const;
        template<> void SmartOption<std::vector<std::string>>::set(std::vector<std::string> str) const;
        template<> void SmartOption<std::vector<double>>::set(std::vector<std::string> str) const;
        template<> void SmartOption<std::vector<int>>::set(std::vector<std::string> str) const;
    }
}