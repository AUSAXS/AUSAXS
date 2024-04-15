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
#include <rigidbody/sequencer/setup/LoadElement.h>
#include <rigidbody/sequencer/setup/ConstraintElement.h>
#include <rigidbody/sequencer/setup/AutoConstraintsElement.h>
#include <rigidbody/constraints/DistanceConstraint.h>

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
    LoadElement,
    Constraint,
    AutomaticConstraint,
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
        {ElementType::LoadElement, {"load", "open"}},
        {ElementType::Constraint, {"constraint", "constrain"}},
        {ElementType::AutomaticConstraint, {"generate_constraints", "autoconstraints", "autoconstrain"}},
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

settings::rigidbody::ConstraintGenerationStrategyChoice get_constraint_strategy(std::string_view line) {
    if (line == "none") {return settings::rigidbody::ConstraintGenerationStrategyChoice::None;}
    if (line == "linear") {return settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;}
    if (line == "volumetric") {return settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric;}
    throw except::invalid_argument("SequenceParser::get_constraint_strategy: Unknown strategy \"" + std::string(line) + "\"");
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

    std::string body1, body2, type;
    int iatom1, iatom2;

    bool found_body1 = false, found_body2 = false, found_iatom1 = false, found_iatom2 = false, found_type = false;
    for (const auto& name : valid_args[Args::body1]) {
        if (args.contains(name)) {
            found_body1 = true;
            body1 = args.at(name)[0];
            break;
        }
    }

    for (const auto& name : valid_args[Args::body2]) {
        if (args.contains(name)) {
            found_body2 = true;
            body2 = args.at(name)[0];
            break;
        }
    }

    for (const auto& name : valid_args[Args::type]) {
        if (args.contains(name)) {
            found_type = true;
            type = args.at(name)[0];
            break;
        }
    }

    for (const auto& name : valid_args[Args::iatom1]) {
        if (args.contains(name)) {
            found_iatom1 = true;
            iatom1 = std::stoi(args.at(name)[0]);
            break;
        }
    }

    for (const auto& name : valid_args[Args::iatom2]) {
        if (args.contains(name)) {
            found_iatom2 = true;
            iatom2 = std::stoi(args.at(name)[0]);
            break;
        }
    }

    if (!found_body1) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"body1\".");}
    if (!found_body2) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"body2\".");}

    if (!(found_iatom1 && found_iatom2)) {
        if (!found_type) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"iatom1\", and \"iatom2\", or \"type\".");}
        if (type == "closest") {
            return std::make_unique<ConstraintElement>(
                static_cast<Sequencer*>(loop_stack.front()), 
                body1, 
                body2
            );
        }
        if (type == "center_mass") {
            return std::make_unique<ConstraintElement>(
                static_cast<Sequencer*>(loop_stack.front()), 
                body1, 
                body2,
                true
            );
        }
        throw except::invalid_argument("SequenceParser::parse_arguments: Invalid argument for \"type\": \"" + type + "\".");
    } else {
        return std::make_unique<ConstraintElement>(
            static_cast<Sequencer*>(loop_stack.front()), 
            body1, 
            body2,
            iatom1,
            iatom2
        );
    }
    throw except::invalid_argument("SequenceParser::parse_arguments: Invalid arguments for \"constraint\".");
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::AutomaticConstraint>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"auto_constraints\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<AutoConstraintsElement>(static_cast<Sequencer*>(loop_stack.front()), get_constraint_strategy(args.begin()->second[0]));
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::LoadElement>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    enum class Args {paths, splits, names};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::paths, {"path", "paths", "load", "anonymous"}},
        {Args::splits, {"splits", "split"}},
        {Args::names, {"names", "name"}}
    };
    if (args.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"load\". Expected at least one argument.");}

    std::vector<std::string> paths, names;
    std::vector<int> splits;

    std::string current_arg;
    for (const auto& name : valid_args[Args::paths]) {
        if (args.contains(name)) {
            current_arg = name;
            paths = args.at(name);
            break;
        }
    }
    if (current_arg.empty()) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"paths\".");}

    for (const auto& name : valid_args[Args::splits]) {
        if (args.contains(name)) {
            for (const auto& split : args.at(name)) {
                try {
                    splits.push_back(std::stoi(split));
                } catch (std::invalid_argument& e) {
                    throw except::invalid_argument("SequenceParser::parse_arguments: Invalid argument for splits: \"" + split + "\".");
                }
            }
            break;
        }
    }

    for (const auto& name : valid_args[Args::names]) {
        if (args.contains(name)) {
            names = args.at(name);
            break;
        }
    }

    if (!splits.empty()) {
        if (paths.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Splits can only be used with a single path.");}
        return std::make_unique<LoadElement>(static_cast<Sequencer*>(loop_stack.front()), paths[0], splits, names);
    }
    return std::make_unique<LoadElement>(static_cast<Sequencer*>(loop_stack.front()), paths, names);
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

    int iterations;
    double radians, angstroms;
    settings::rigidbody::ParameterGenerationStrategyChoice strategy;
    settings::rigidbody::DecayStrategyChoice decay_strategy;

    std::string current_arg;
    try {
        for (const auto& name : valid_args[Args::iterations]) {
            if (args.contains(name)) {
                current_arg = name;
                iterations = std::stoi(args.at(name)[0]);
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
                angstroms = std::stoi(args.at(name)[0]);
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
                radians = std::stod(args.at(name)[0]);
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
            rigidbody::factory::create_decay_strategy(iterations, decay_strategy),
            angstroms,
            radians,
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
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Save>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"save\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<SaveElement>(loop_stack.back(), args.begin()->second[0]);
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::SaveOnImprove>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"save_on_improve\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    auto last_element_cast = dynamic_cast<OptimizeStepElement*>(loop_stack.back()->_get_elements().back().get());
    if (last_element_cast == nullptr) {
        throw except::invalid_argument("SequenceParser::parse_arguments: \"save_on_improve\" must be called after an \"optimize_step\" element.");
    }
    last_element_cast->save_on_improvement(args.begin()->second[0]);
    return nullptr;
}

std::unique_ptr<Sequencer> SequenceParser::parse(const io::ExistingFile& config, const io::ExistingFile& saxs) {
    std::ifstream in(config.path());
    if (!in.is_open()) {throw except::io_error("SequenceParser::parse: Could not open file \"" + config.path() + "\".");}

    std::unique_ptr<Sequencer> sequencer = std::make_unique<Sequencer>(saxs);
    loop_stack = {sequencer.get()};

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
                for (int i = 0; i < sub_tokens.size(); i++) {
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

            case ElementType::LoadElement:
                loop_stack.back()->_get_elements().push_back(parse_arguments<ElementType::LoadElement>(args));
                break;

            default:
                throw except::invalid_argument("SequenceParser::constructor: Unknown element type.");
        }
    }
    return sequencer;
}