// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/sequencer/TransformElement.h>
#include <rigidbody/sequencer/OptimizeStepElement.h>
#include <rigidbody/sequencer/EveryNStepElement.h>
#include <rigidbody/sequencer/OnImprovementElement.h>
#include <rigidbody/sequencer/SaveElement.h>
#include <rigidbody/sequencer/setup/LoadElement.h>
#include <rigidbody/sequencer/setup/ConstraintElement.h>
#include <rigidbody/sequencer/setup/AutoConstraintsElement.h>
#include <rigidbody/sequencer/setup/RelativeHydrationElement.h>
#include <rigidbody/sequencer/setup/SymmetryElement.h>

#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/BodySplitter.h>

#include <utility/observer_ptr.h>
#include <utility/StringUtils.h>
#include <utility/Exceptions.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <io/ExistingFile.h>
#include <settings/RigidBodySettings.h>

#include <fstream>
#include <unordered_map>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

enum class rigidbody::sequencer::ElementType {
    OverlapStrength,
    LoadElement,
    SymmetryElement,
    Constraint,
    AutomaticConstraint,
    LoopBegin,
    LoopEnd,
    Parameter,
    BodySelect,
    Transform,
    OptimizeStep,
    EveryNStep,
    OnImprovement,
    Save,
    RelativeHydration
};

ElementType get_type(std::string_view line) {
    static std::unordered_map<ElementType, std::vector<std::string>> type_map = {
        {ElementType::OverlapStrength, {"overlap_strength"}},
        {ElementType::LoadElement, {"load", "open"}},
        {ElementType::SymmetryElement, {"symmetry"}},
        {ElementType::Constraint, {"constraint", "constrain"}},
        {ElementType::AutomaticConstraint, {"generate_constraints", "autoconstraints", "autoconstrain"}},
        {ElementType::LoopBegin, {"loop"}},
        {ElementType::LoopEnd, {"end"}},
        {ElementType::Parameter, {"parameter"}},
        {ElementType::BodySelect, {"selection", "body_selection"}},
        {ElementType::Transform, {"transform"}},
        {ElementType::OptimizeStep, {"optimize_step", "optimize_once"}},
        {ElementType::EveryNStep, {"every", "for_every"}},
        {ElementType::OnImprovement, {"on_improvement"}},
        {ElementType::Save, {"save", "write"}},
        {ElementType::RelativeHydration, {"relative_hydration"}}
    };
    for (const auto& [type, prefixes] : type_map) {
        for (const auto& prefix : prefixes) {
            if (line.starts_with(prefix)) {return type;}
        }
    }
    throw except::invalid_argument("SequenceParser::get_type: Unknown element type \"" + std::string(line) + "\".");
}

settings::rigidbody::TransformationStrategyChoice get_transform_strategy(std::string_view line) {
    if (line == "rigid_transform") {return settings::rigidbody::TransformationStrategyChoice::RigidTransform;}
    if (line == "single_transform") {return settings::rigidbody::TransformationStrategyChoice::SingleTransform;}
    throw except::invalid_argument("SequenceParser::get_transform_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

settings::rigidbody::BodySelectStrategyChoice get_body_select_strategy(std::string_view line) {
    if (line == "random_body") {return settings::rigidbody::BodySelectStrategyChoice::RandomBodySelect;}
    if (line == "random_constraint") {return settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect;}
    if (line == "sequential_body") {return settings::rigidbody::BodySelectStrategyChoice::SequentialBodySelect;}
    if (line == "sequential_constraint") {return settings::rigidbody::BodySelectStrategyChoice::SequentialConstraintSelect;}
    throw except::invalid_argument("SequenceParser::get_body_select_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

settings::rigidbody::DecayStrategyChoice get_decay_strategy(std::string_view line) {
    if (line == "linear") {return settings::rigidbody::DecayStrategyChoice::Linear;}
    if (line == "exponential") {return settings::rigidbody::DecayStrategyChoice::Exponential;}
    if (line == "none") {return settings::rigidbody::DecayStrategyChoice::None;}
    throw except::invalid_argument("SequenceParser::get_decay_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

settings::rigidbody::ParameterGenerationStrategyChoice get_parameter_strategy(std::string_view line) {
    if (line == "rotate_only") {return settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly;}
    if (line == "translate_only") {return settings::rigidbody::ParameterGenerationStrategyChoice::TranslationsOnly;}
    if (line == "symmetry_only") {return settings::rigidbody::ParameterGenerationStrategyChoice::SymmetryOnly;}
    if (line == "both" || line == "rotate_and_translate") {return settings::rigidbody::ParameterGenerationStrategyChoice::Simple;}
    throw except::invalid_argument("SequenceParser::get_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

settings::rigidbody::ConstraintGenerationStrategyChoice get_constraint_strategy(std::string_view line) {
    if (line == "none") {return settings::rigidbody::ConstraintGenerationStrategyChoice::None;}
    if (line == "linear") {return settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;}
    if (line == "volumetric") {return settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric;}
    throw except::invalid_argument("SequenceParser::get_constraint_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

template<typename T>
struct ArgResult {T value; bool found;};

template<typename T>
ArgResult<T> check_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args);

template<>
ArgResult<std::string> check_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            return {args.at(name)[0], true};
        }
    }
    return {"", false};
}

template<>
ArgResult<std::vector<std::string>> check_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            return {args.at(name), true};
        }
    }
    return {{}, false};
}

template<>
ArgResult<int> check_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            try {
                return {std::stoi(args.at(name)[0]), true};
            } catch (std::exception&) {
                throw except::invalid_argument("SequenceParser::check_arg: \"" + args.at(name)[0] + "\" cannot be interpreted as an integer.");
            }
        }
    }
    return {0, false};
}

template<>
ArgResult<std::vector<int>> check_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            std::vector<int> values;
            try {
                for (const auto& value : args.at(name)) {values.push_back(std::stoi(value));}
            } catch (std::exception&) {
                throw except::invalid_argument("SequenceParser::check_arg: \"" + args.at(name)[values.size()] + "\" cannot be interpreted as an integer.");
            }
            return {std::move(values), true};
        }
    }
    return {{}, false};
}

template<>
ArgResult<double> check_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            try {
                return {std::stod(args.at(name)[0]), true};
            } catch (std::exception&) {
                throw except::invalid_argument("SequenceParser::check_arg: \"" + args.at(name)[0] + "\" cannot be interpreted as a decimal value.");
            }
        }
    }
    return {0.0, false};
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Constraint>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    enum class Args {body1, body2, iatom1, iatom2, type};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::body1, {"first", "first_body"}},
        {Args::body2, {"second", "second_body"}},
        {Args::iatom1, {"iatom1", "iatom_1", "i1"}},
        {Args::iatom2, {"iatom2", "iatom_2", "i2"}},
        {Args::type, {"type", "kind"}}
    };
    if (args.size() < 3) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"constraint\". Expected at least 3, but got " + std::to_string(args.size()) + ".");}

    auto body1 = check_arg<std::string>(valid_args[Args::body1], args);
    auto body2 = check_arg<std::string>(valid_args[Args::body2], args);
    auto type  = check_arg<std::string>(valid_args[Args::type], args);
    auto iatom1 = check_arg<int>(valid_args[Args::iatom1], args);
    auto iatom2 = check_arg<int>(valid_args[Args::iatom2], args);

    if (!body1.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"body1\".");}
    if (!body1.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"body2\".");}

    if (!(iatom1.found && iatom2.found)) {
        if (!type.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"iatom1\", and \"iatom2\", or \"type\".");}
        if (type.value == "closest") {
            return std::make_unique<ConstraintElement>(
                static_cast<Sequencer*>(loop_stack.front()), 
                body1.value, 
                body2.value
            );
        }
        if (type.value == "center_mass") {
            return std::make_unique<ConstraintElement>(
                static_cast<Sequencer*>(loop_stack.front()), 
                body1.value, 
                body2.value,
                true
            );
        }
        throw except::invalid_argument("SequenceParser::parse_arguments: Invalid argument for \"type\": \"" + type.value + "\".");
    } else {
        return std::make_unique<ConstraintElement>(
            static_cast<Sequencer*>(loop_stack.front()), 
            body1.value, 
            body2.value,
            iatom1.value,
            iatom2.value
        );
    }
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::AutomaticConstraint>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"auto_constraints\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<AutoConstraintsElement>(static_cast<Sequencer*>(loop_stack.front()), get_constraint_strategy(args.begin()->second[0]));
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::LoadElement>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    enum class Args {paths, splits, names, saxs};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::paths, {"path", "paths", "load", "pdb", "anonymous"}},
        {Args::splits, {"splits", "split"}},
        {Args::names, {"names", "name"}},
        {Args::saxs, {"saxs_path", "saxs"}}
    };
    if (args.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: \"load\": Invalid number of arguments. Expected at least one.");}

    auto paths = check_arg<std::vector<std::string>>(valid_args[Args::paths], args);
    auto names = check_arg<std::vector<std::string>>(valid_args[Args::names], args);
    auto saxs_path = check_arg<std::string>(valid_args[Args::saxs], args);
    if (!paths.found) {throw except::invalid_argument("SequenceParser::parse_arguments: \"load\": Missing required argument \"paths\".");}

    {   // check the special case of splitting by chainID
        auto split_test = check_arg<std::string>(valid_args[Args::splits], args);
        if (split_test.found && split_test.value == "chain") {
            if (paths.value.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: \"load\": Chain splitting can only be used with a single path.");}
            return std::make_unique<LoadElement>(static_cast<Sequencer*>(loop_stack.front()), paths.value[0], names.value, saxs_path.value);
        }
    }

    auto splits = check_arg<std::vector<int>>(valid_args[Args::splits], args);
    if (splits.found) {
        if (paths.value.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: \"load\": Splits can only be used with a single path.");}
        return std::make_unique<LoadElement>(static_cast<Sequencer*>(loop_stack.front()), paths.value[0], splits.value, names.value, saxs_path.value);
    }
    return std::make_unique<LoadElement>(static_cast<Sequencer*>(loop_stack.front()), paths.value, names.value, saxs_path.value);
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::SymmetryElement>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    auto rigidbody = loop_stack.front()->_get_rigidbody();
    if (rigidbody->size_body() < args.size()) {
        throw except::invalid_argument(
            "SequenceParser::parse_arguments: Too many arguments for \"symmetry\". "
            "Expected no more than " + std::to_string(rigidbody->size_body()) + "."
        );
    }

    // anonymous arg support
    if (args.size() == 1) {
        return std::make_unique<SymmetryElement>(static_cast<Sequencer*>(loop_stack.front()), std::vector<std::string>{"b1"}, std::vector{symmetry::get(args.begin()->second[0])});
    }

    std::vector<symmetry::type> symmetries;
    std::vector<std::string> names;
    for (const auto& [name, value] : args) {
        names.push_back(name);
        symmetries.push_back(symmetry::get(value[0]));
    }
    return std::make_unique<SymmetryElement>(static_cast<Sequencer*>(loop_stack.front()), names, symmetries);
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::LoopBegin>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"loop\". Expected 1, but got " + std::to_string(args.size()) + ".");}

    int iterations;
    try {
        iterations = std::stoi(args.begin()->second[0]);
    } catch (std::invalid_argument& e) {
        throw except::invalid_argument("SequenceParser::parse_arguments: Invalid argument for loop iterations: \"" + args.begin()->second[0] + "\".");
    }
    return std::make_unique<LoopElement>(
        loop_stack.back(), 
        iterations
    );
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Parameter>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    enum class Args {iterations, angstroms, radians, strategy, decay_strategy};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::iterations, {"iterations", "decay_over"}},
        {Args::angstroms, {"angstroms", "translation_amplitude", "translate"}},
        {Args::radians, {"radians", "rotation_amplitude", "rotate"}},
        {Args::strategy, {"strategy"}},
        {Args::decay_strategy, {"decay_strategy", "decay"}}
    };
    if (args.size() != valid_args.size()) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"parameter\". Expected 5, but got " + std::to_string(args.size()) + ".");}

    settings::rigidbody::ParameterGenerationStrategyChoice strategy = settings::rigidbody::ParameterGenerationStrategyChoice::Simple;
    settings::rigidbody::DecayStrategyChoice decay_strategy = settings::rigidbody::DecayStrategyChoice::Linear;

    auto iterations = check_arg<int>(valid_args[Args::iterations], args);
    auto angstroms = check_arg<double>(valid_args[Args::angstroms], args);
    auto radians = check_arg<double>(valid_args[Args::radians], args);

    if (!iterations.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"iterations\".");}
    if (!angstroms.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"angstroms\".");}
    if (!radians.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"radians\".");}

    std::string current_arg = "";
    for (const auto& name : valid_args[Args::strategy]) {
        if (args.contains(name)) {
            current_arg = name;
            strategy = get_parameter_strategy(args.at(name)[0]);
            break;
        }
    }
    if (current_arg.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"strategy\".");}

    current_arg = "";
    for (const auto& name : valid_args[Args::decay_strategy]) {
        if (args.contains(name)) {
            current_arg = name;
            decay_strategy = get_decay_strategy(args.at(name)[0]);
            break;
        }
    }
    if (current_arg.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"decay_strategy\".");}

    return std::make_unique<ParameterElement>(
        loop_stack.back(),
        rigidbody::factory::create_parameter_strategy(
            loop_stack.front()->_get_rigidbody(),
            rigidbody::factory::create_decay_strategy(iterations.value, decay_strategy),
            angstroms.value,
            radians.value,
            strategy
        )
    );
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::BodySelect>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"body_select\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    settings::rigidbody::BodySelectStrategyChoice strategy = get_body_select_strategy(args.begin()->second[0]);
    return std::make_unique<BodySelectElement>(
        loop_stack.back(),
        rigidbody::factory::create_selection_strategy(
            loop_stack.front()->_get_rigidbody(),
            strategy
        )
    );
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::LoopEnd>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 0) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"end\". Expected 0, but got " + std::to_string(args.size()) + ".");}
    return nullptr;
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Transform>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"transform\". Expected 1, but got " + std::to_string(args.size()) + ".");}

    settings::rigidbody::TransformationStrategyChoice strategy = get_transform_strategy(args.begin()->second[0]);
    return std::make_unique<TransformElement>(
        loop_stack.back(),
        rigidbody::factory::create_transform_strategy(
            loop_stack.front()->_get_rigidbody(),
            strategy
        )
    );
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::RelativeHydration>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    enum class Options {Maximum, High, Normal, Low, Minimum};
    std::unordered_map<std::string, Options> options = {
        {"max",     Options::Maximum},
        {"maximum", Options::Maximum},
        {"high",    Options::High},
        {"normal",  Options::Normal},
        {"low",     Options::Low},
        {"minimum", Options::Minimum},
        {"min",     Options::Minimum}
    };

    auto to_value = [] (Options opt) {
        switch (opt) {
            case Options::Maximum:  return 1.75;
            case Options::High:     return 1.5;
            case Options::Normal:   return 1.0;
            case Options::Low:      return 0.5;
            case Options::Minimum:  return 0.25;
        }
        return 0.0;
    };

    auto rigidbody = loop_stack.front()->_get_rigidbody();
    if (args.size() != rigidbody->size_body()) {
        throw except::invalid_argument(
            "SequenceParser::parse_arguments: Invalid number of arguments for \"relative_hydration\". "
            "Expected " + std::to_string(rigidbody->size_body()) + ", but got " + std::to_string(args.size()) + "."
        );
    }

    std::vector<double> ratios;
    std::vector<std::string> names;
    for (const auto& [name, value] : args) {
        names.push_back(name);
        ratios.push_back(to_value(options[value[0]]));
    }
    return std::make_unique<RelativeHydrationElement>(static_cast<Sequencer*>(loop_stack.front()), names, ratios);
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::OptimizeStep>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 0) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"optimize_step\". Expected 0, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<OptimizeStepElement>(loop_stack.back());
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::EveryNStep>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"every_n_step\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<EveryNStepElement>(loop_stack.back(), std::stoi(args.begin()->second[0]));
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::OnImprovement>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (!args.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"on_improvement\". Expected 0, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<OnImprovementElement>(loop_stack.back());
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Save>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"save\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<SaveElement>(loop_stack.back(), args.begin()->second[0]);
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::OverlapStrength>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    enum class Args {scaling, distance};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::scaling, {"scaling", "factor"}},
        {Args::distance, {"max", "max_distance", "distance"}}
    };
    if (args.size() != valid_args.size()) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"overlap_strength\". Expected 2, but got " + std::to_string(args.size()) + ".");}

    auto scaling  = check_arg<double>(valid_args[Args::scaling], args);
    auto distance = check_arg<double>(valid_args[Args::distance], args);
    if (!scaling.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"scaling\".");}
    if (!distance.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"distance\".");}

    static_cast<Sequencer*>(loop_stack.front())->set_overlap_function([a=scaling.value, d=distance.value] (double x) {return x < d ? a*std::pow((d-x)/d, 2) : 0;});
    return nullptr;
}

std::unique_ptr<Sequencer> SequenceParser::parse(const io::ExistingFile& config) {
    std::ifstream in(config.path());
    if (!in.is_open()) {throw except::io_error("SequenceParser::parse: Could not open file \"" + config.path() + "\".");}

    // the main sequencer object
    std::unique_ptr<Sequencer> sequencer = std::make_unique<Sequencer>();

    // the loop stack
    // the top element of this stack is the current loop element which new elements will be added to
    // note that the sequencer itself is just a dummy loop element with an iteration count of 1
    loop_stack = {sequencer.get()};
    sequencer->_set_config_folder(config.directory());
    
    std::string line;
    while(!in.eof()) {
        // read the next line
        std::getline(in, line);
        if (line.empty()) {continue;}

        // skip empty lines & comments
        std::unordered_map<std::string, std::vector<std::string>> args;
        auto tokens = utility::split(line, " \t\r\n");
        if (tokens.empty() || tokens.front()[0] == '#') {continue;}

        // check if the argument list spans over multiple lines
        if (line.find_first_of("{[(") != std::string::npos) {
            // parse args in the next lines
            std::string argline;
            std::getline(in, argline);
            while (argline.find_first_of("]})") == std::string::npos) {
                if (argline.empty()) {continue;}
                auto sub_tokens = utility::split(argline, " \t\r\n");

                // check for comment tokens 
                int end = sub_tokens.size();
                for (unsigned int i = 0; i < sub_tokens.size(); i++) {
                    if (sub_tokens[i][0] == '#') {
                        end = i;
                        break;
                    }
                }

                args[sub_tokens[0]] = std::vector<std::string>(sub_tokens.begin()+1, sub_tokens.begin()+end);
                std::getline(in, argline);
                if (in.eof()) {throw except::io_error("SequenceParser::parse: Unescaped argument list starting in line \"" + line + "\".");}
            }
        } 
        
        // else check if we only have a single argument
        else if (tokens.size() == 2) {
            // allow for a single anonymous argument
            args["anonymous"] = {tokens[1]};
        }

        std::cout << "continuing with element: " << std::endl;
        for (const auto& t : tokens) {
            std::cout << "\t\"" << t << "\"" << std::endl;
        }
        std::cout << "and arguments: " << std::endl; 
        for (const auto& [key, value] : args) {
            std::cout << "\t\"" << key << "\"" << std::endl;
            for (const auto& v : value) {
                std::cout << "\t\t\"" << v << "\"" << std::endl;
            }
        }
        switch (get_type(tokens[0])) {
            case ElementType::Constraint:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::Constraint>(args));
                break;

            case ElementType::AutomaticConstraint:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::AutomaticConstraint>(args));
                break;

            case ElementType::LoopBegin:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::LoopBegin>(args));
                loop_stack.push_back(static_cast<LoopElement*>(loop_stack.back()->_get_elements().back().get()));
                break;

            case ElementType::LoopEnd:
                parse_arguments<ElementType::LoopEnd>(args);
                if (loop_stack.size() == 1) {throw except::invalid_argument("SequenceParser::constructor: Unmatched \"end\" statement.");}
                loop_stack.pop_back();
                break;

            case ElementType::Parameter:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::Parameter>(args));
                break;

            case ElementType::BodySelect:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::BodySelect>(args));
                break;

            case ElementType::Transform:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::Transform>(args));
                break;
            
            case ElementType::OptimizeStep:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::OptimizeStep>(args));
                break;

            case ElementType::EveryNStep:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::EveryNStep>(args));
                loop_stack.push_back(static_cast<LoopElement*>(loop_stack.back()->_get_elements().back().get()));
                break;

            case ElementType::OnImprovement:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::OnImprovement>(args));
                loop_stack.push_back(static_cast<LoopElement*>(loop_stack.back()->_get_elements().back().get()));
                break;

            case ElementType::Save:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::Save>(args));
                break;

            case ElementType::LoadElement:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::LoadElement>(args));
                break;

            case ElementType::SymmetryElement:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::SymmetryElement>(args));
                break;

            case ElementType::OverlapStrength:
                parse_arguments<ElementType::OverlapStrength>(args);
                break;

            case ElementType::RelativeHydration:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::RelativeHydration>(args));
                break;

            default:
                throw except::invalid_argument("SequenceParser::constructor: Unknown element type \"" + tokens[0] + "\".");
        }
    }
    return sequencer;
}