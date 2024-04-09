/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/sequencer/TransformElement.h>
#include <rigidbody/sequencer/OptimizeStepElement.h>
#include <rigidbody/sequencer/EveryNStepElement.h>
#include <rigidbody/sequencer/SaveElement.h>

#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/transform/TransformFactory.h>

#include <rigidbody/RigidBody.h>
#include <data/BodySplitter.h>
#include <utility/observer_ptr.h>
#include <utility/StringUtils.h>
#include <utility/Exceptions.h>
#include <io/ExistingFile.h>
#include <settings/RigidBodySettings.h>

#include <fstream>
#include <unordered_map>

using namespace rigidbody::sequencer;

enum class rigidbody::sequencer::ElementType {
    LoopBegin,
    LoopEnd,
    Parameter,
    BodySelect,
    Transform,
    OptimizeStep,
    EveryNStep,
    Save,
    SaveOnImprove
};

ElementType get_type(std::string_view line) {
    static std::unordered_map<ElementType, std::vector<std::string>> type_map = {
        {ElementType::LoopBegin, {"loop"}},
        {ElementType::LoopEnd, {"end"}},
        {ElementType::Parameter, {"parameter"}},
        {ElementType::BodySelect, {"selection", "body_selection"}},
        {ElementType::Transform, {"transform"}},
        {ElementType::OptimizeStep, {"optimize_step", "optimize_once"}},
        {ElementType::EveryNStep, {"every", "for_every"}},
        {ElementType::Save, {"save", "write"}},
        {ElementType::SaveOnImprove, {"save_on_improve", "write_on_improve"}}
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
    if (line == "random_body") {return settings::rigidbody::BodySelectStrategyChoice::RandomSelect;}
    if (line == "random_constraint") {return settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect;}
    if (line == "sequential") {return settings::rigidbody::BodySelectStrategyChoice::SequentialSelect;}
    throw except::invalid_argument("SequenceParser::get_body_select_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

settings::rigidbody::DecayStrategyChoice get_decay_strategy(std::string_view line) {
    if (line == "linear") {return settings::rigidbody::DecayStrategyChoice::Linear;}
    if (line == "exponential") {return settings::rigidbody::DecayStrategyChoice::Exponential;}
    throw except::invalid_argument("SequenceParser::get_decay_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

settings::rigidbody::ParameterGenerationStrategyChoice get_parameter_strategy(std::string_view line) {
    if (line == "rotate_only") {return settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly;}
    if (line == "translate_only") {return settings::rigidbody::ParameterGenerationStrategyChoice::TranslationsOnly;}
    if (line == "both" || line == "rotate_and_translate") {return settings::rigidbody::ParameterGenerationStrategyChoice::Simple;}
    throw except::invalid_argument("SequenceParser::get_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::LoopBegin>(const std::unordered_map<std::string, std::string>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"loop\". Expected 1, but got " + std::to_string(args.size()) + ".");}

    int iterations;
    try {
        iterations = std::stoi(args.begin()->second);
    } catch (std::invalid_argument& e) {
        throw except::invalid_argument("SequenceParser::parse_arguments: Invalid argument for loop iterations: \"" + args.begin()->second + "\".");
    }
    return std::make_unique<LoopElement>(
        loop_stack.back(), 
        iterations
    );
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Parameter>(const std::unordered_map<std::string, std::string>& args) {
    enum class Args {iterations, angstroms, radians, strategy, decay_strategy};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::iterations, {"iterations", "decay_over"}},
        {Args::angstroms, {"angstroms", "translation_amplitude", "translate"}},
        {Args::radians, {"radians", "rotation_amplitude", "rotate"}},
        {Args::strategy, {"strategy"}},
        {Args::decay_strategy, {"decay_strategy", "decay"}}
    };
    if (args.size() != valid_args.size()) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"parameter\". Expected 5, but got " + std::to_string(args.size()) + ".");}

    int iterations;
    double radians, angstroms;
    settings::rigidbody::ParameterGenerationStrategyChoice strategy;
    settings::rigidbody::DecayStrategyChoice decay_strategy;

    std::string current_arg;
    try {
        for (const auto& name : valid_args[Args::iterations]) {
            if (args.contains(name)) {
                current_arg = name;
                iterations = std::stoi(args.at(name));
                break;
            }
        }
    } catch (std::invalid_argument& e) {
        throw except::invalid_argument("SequenceParser::parse_arguments: Invalid argument for parameter iterations: \"" + current_arg + "\".");
    }
    if (current_arg.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"iterations\".");}

    current_arg = "";
    try {
        for (const auto& name : valid_args[Args::angstroms]) {
            if (args.contains(name)) {
                current_arg = name;
                angstroms = std::stoi(args.at(name));
                break;
            }
        }
    } catch (std::invalid_argument& e) {
        throw except::invalid_argument("SequenceParser::parse_arguments: Invalid argument for parameter angstroms: \"" + current_arg + "\".");
    }
    if (current_arg.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"angstroms\".");}

    current_arg = "";
    try {
        for (const auto& name : valid_args[Args::radians]) {
            if (args.contains(name)) {
                current_arg = name;
                radians = std::stod(args.at(name));
                break;
            }
        }
    } catch (std::invalid_argument& e) {
        throw except::invalid_argument("SequenceParser::parse_arguments: Invalid argument for parameter radians: \"" + current_arg + "\".");
    }
    if (current_arg.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"radians\".");}

    current_arg = "";
    for (const auto& name : valid_args[Args::strategy]) {
        if (args.contains(name)) {
            current_arg = name;
            strategy = get_parameter_strategy(args.at(name));
            break;
        }
    }
    if (current_arg.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"strategy\".");}

    current_arg = "";
    for (const auto& name : valid_args[Args::decay_strategy]) {
        if (args.contains(name)) {
            current_arg = name;
            decay_strategy = get_decay_strategy(args.at(name));
            break;
        }
    }
    if (current_arg.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"decay_strategy\".");}

    return std::make_unique<ParameterElement>(
        loop_stack.back(),
        rigidbody::factory::create_parameter_strategy(
            rigidbody::factory::create_decay_strategy(iterations, decay_strategy),
            angstroms,
            radians,
            strategy
        )
    );
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::BodySelect>(const std::unordered_map<std::string, std::string>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"body_select\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    settings::rigidbody::BodySelectStrategyChoice strategy = get_body_select_strategy(args.begin()->second);
    return std::make_unique<BodySelectElement>(
        loop_stack.back(),
        rigidbody::factory::create_selection_strategy(
            loop_stack.front()->_get_rigidbody(),
            strategy
        )
    );
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::LoopEnd>(const std::unordered_map<std::string, std::string>& args) {
    if (args.size() != 0) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"end\". Expected 0, but got " + std::to_string(args.size()) + ".");}
    return nullptr;
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Transform>(const std::unordered_map<std::string, std::string>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"transform\". Expected 1, but got " + std::to_string(args.size()) + ".");}

    settings::rigidbody::TransformationStrategyChoice strategy = get_transform_strategy(args.begin()->second);
    return std::make_unique<TransformElement>(
        loop_stack.back(),
        rigidbody::factory::create_transform_strategy(
            loop_stack.front()->_get_rigidbody(),
            strategy
        )
    );
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::OptimizeStep>(const std::unordered_map<std::string, std::string>& args) {
    if (args.size() != 0) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"optimize_step\". Expected 0, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<OptimizeStepElement>(loop_stack.back());
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::EveryNStep>(const std::unordered_map<std::string, std::string>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"every_n_step\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<EveryNStepElement>(loop_stack.back(), std::stoi(args.begin()->second));
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Save>(const std::unordered_map<std::string, std::string>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"save\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<SaveElement>(loop_stack.back(), args.begin()->second);
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::SaveOnImprove>(const std::unordered_map<std::string, std::string>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"save_on_improve\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    auto last_element_cast = dynamic_cast<OptimizeStepElement*>(loop_stack.back()->_get_elements().back().get());
    if (last_element_cast == nullptr) {
        throw except::invalid_argument("SequenceParser::parse_arguments: \"save_on_improve\" must be called after an \"optimize_step\" element.");
    }
    last_element_cast->save_on_improvement(args.begin()->second);
    return nullptr;
}

std::unique_ptr<Sequencer> SequenceParser::parse(const io::ExistingFile& config, const io::ExistingFile& saxs, observer_ptr<RigidBody> rigidbody) {
    std::ifstream in(config.path());
    if (!in.is_open()) {throw except::io_error("SequenceParser::parse: Could not open file \"" + config.path() + "\".");}

    std::unique_ptr<Sequencer> sequencer = std::make_unique<Sequencer>(saxs, rigidbody);
    loop_stack = {sequencer.get()};

    std::string line;
    while(!in.eof()) {
        // verify that the line is a valid element
        std::getline(in, line);
        if (line.empty()) {continue;}
        std::cout << line << std::endl;

        std::unordered_map<std::string, std::string> args;
        auto tokens = utility::split(line, " \t");
        if (line.find_first_of("{[(") != std::string::npos) {
            // parse args in the next lines
            std::string argline;
            std::getline(in, argline);
            while (argline.find_first_of("]})") == std::string::npos) {
                if (argline.empty()) {continue;}
                auto sub_tokens = utility::split(argline, " \t");
                if (sub_tokens.size() != 2) {throw except::invalid_argument("SequenceParser::parse: Invalid argument line \"" + argline + "\".");}
                args[sub_tokens[0]] = sub_tokens[1];
                std::getline(in, argline);
                if (in.eof()) {throw except::io_error("SequenceParser::parse: Unescaped argument list starting in line \"" + line + "\".");}
            }
        } else if (tokens.size() == 2) {
            // allow for a single anonymous argument
            args["anonymous"] = tokens[1];
        }

        for (const auto& [key, value] : args) {
            std::cout << "\t" << key << " = " << value << std::endl;
        }
        switch (get_type(tokens[0])) {
            case ElementType::LoopBegin:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::LoopBegin>(args));
                loop_stack.push_back(static_cast<LoopElement*>(loop_stack.back()->_get_elements().back().get()));
                break;

            case ElementType::LoopEnd:
                parse_arguments<ElementType::LoopEnd>(args);
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
                break;

            case ElementType::Save:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::Save>(args));
                break;

            case ElementType::SaveOnImprove:
                parse_arguments<ElementType::SaveOnImprove>(args);
                break;

            default:
                throw except::invalid_argument("SequenceParser::constructor: Unknown element type.");
        }
    }
    return sequencer;
}