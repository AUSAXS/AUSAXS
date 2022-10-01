#pragma once

#include <vector>
#include <string>
#include <memory>

#include <math/Vector3.h>
#include <utility/Axis.h>
#include <utility/Type.h>
#include <utility/Exceptions.h>

// A small container of the various settings. These should be set *before* their respective classes are instantiated. 
namespace setting {
    struct general {
        static constexpr char residue_folder[] = "temp/residues/"; // Download location for all ligand files. Must be constexpr. 
    };

    struct figures {
        inline static std::string format = "pdf"; // The output format.
    };

    struct grid {
        enum PlacementStrategy {AxesStrategy, RadialStrategy, JanStrategy};
        enum CullingStrategy {CounterStrategy, OutlierStrategy, RandomStrategy};

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
        inline static double scattering_intensity_plot_binned_width = 0.1; // The width of each bin for the scattering plots.
        inline static unsigned int bins = 1000;                            // The number of bins for the scattering plots.
        inline static double qmin = 0.001;                                 // Lower limit on the used q-values
        inline static double qmax = 1.001;                                 // Upper limit on the used q-values
    };

    struct fit {
        inline static unsigned int N = 100; // Number of points sampled when discretizing a model scattering curve
        inline static bool verbose = false; // Whether to print the fit progress to the console
    };

    struct rigidbody {
        enum TransformationStrategyChoice {RigidTransform};
        enum ParameterGenerationStrategyChoice {Simple};
        enum BodySelectStrategyChoice {RandomSelect, SequentialSelect};

        inline static TransformationStrategyChoice tsc = TransformationStrategyChoice::RigidTransform;
        inline static ParameterGenerationStrategyChoice pgsc = ParameterGenerationStrategyChoice::Simple;
        inline static BodySelectStrategyChoice bssc = BodySelectStrategyChoice::RandomSelect;
    };

    struct em {
        enum CullingStrategyChoice {NoStrategy, CounterStrategy};
        CullingStrategyChoice csc = CullingStrategyChoice::CounterStrategy; // The choice of culling algorithm. 

        inline static unsigned int sample_frequency = 1; // How often a bin is sampled in any direction. 
        inline static double concentration = 2;          // The concentration in mg/mL used when calculating the absolute intensity scale for simulations.
        inline static unsigned int charge_levels = 100;  // The number of partial histograms to utilize.

        inline static bool hydrate = true;               // Whether to hydrate the protein in the EM algorithm.
        inline static unsigned int evals = 100;          // Base number of evaluations used in the EM fitter. 

        struct simulation {
            inline static bool noise = true; // Whether to generate noise for the simulations. 
        };
    };

    struct plot {
        inline static std::string path = "figures/"; // The path to the output folder. 

        struct image {
            inline static std::vector<double> contour = {}; // The contour levels for the image plots.
        };

        struct em {
            // Whether to plot the evaluated chi2 points. 
            // Produces 2 plots; one of the full landscape and another of the area near the minimum. 
            // The number of points is roughly determined by setting::em::evals
            inline static bool plot_cutoff_points = true; 
        };
    };

    void read(const std::string path);
    void write(const std::string path);

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
            void set(std::vector<std::string>) const override {throw except::unexpected("Error in \"Settings::SmartOption::set\": Not implemented for type \"" + type(setting) + "\".");}

            /**
             * @brief Get a string representation of the setting represented by this SmartOption.
             *        This method must be specialized for all possible setting types.
             */
            std::string get() const override {throw except::unexpected("Error in \"Settings::SmartOption::set\": Not implemented for type \"" + type(setting) + "\".");}
        };

        /**
         * @brief Helper function for creating shared pointers to SmartOptions. This is to automatically deduce template types,
         *        which is not possible with the default shared_ptr constructor.
         */
        template<typename T>
        std::shared_ptr<SmartOption<T>> make_shared(std::vector<std::string> aliases, T& setting) {
            return std::make_shared<SmartOption<T>>(aliases, setting);
        }

        inline static const std::vector<std::shared_ptr<ISmartOption>> options = {
            // figures
            make_shared({"format"}, setting::figures::format),

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
            make_shared({"scattering-intensity-plot-binned-width"}, setting::axes::scattering_intensity_plot_binned_width),
            make_shared({"qlow", "qmin"}, setting::axes::qmin),
            make_shared({"qhigh", "qmax"}, setting::axes::qmax),

            // fit
            make_shared({"N"}, setting::fit::N),

            // rigidbody

            // em
            make_shared({"sample-frequency"}, setting::em::sample_frequency),
            make_shared({"concentration"}, setting::em::concentration),
            make_shared({"charge-levels"}, setting::em::charge_levels),
            make_shared({"noise"}, setting::em::simulation::noise),
            make_shared({"contour"}, setting::plot::image::contour),

            // plot

        };

        // declare template specializations
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
    }
}