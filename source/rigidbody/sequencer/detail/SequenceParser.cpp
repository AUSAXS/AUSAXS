// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/ConstraintFactory.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/Rigidbody.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <utility/observer_ptr.h>
#include <utility/StringUtils.h>
#include <utility/Exceptions.h>
#include <utility/Logging.h>
#include <utility/Random.h>
#include <io/ExistingFile.h>
#include <settings/RigidBodySettings.h>
#include <settings/GeneralSettings.h>

#include <fstream>
#include <unordered_map>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

enum class rigidbody::sequencer::ElementType {
    AutomaticConstraint,
    BodySelect,
    Constraint,
    Copy,
    EveryNStep,
    LoadElement,
    Log,
    LoopBegin,
    LoopEnd,
    Message,
    OnImprovement,
    OptimizeStep,
    OutputFolder,
    OverlapStrength,
    Parameter,
    RelativeHydration,
    Save,
    Seed,
    SymmetryElement,
    Transform,
};

ElementType get_type(std::string_view line) {
    static std::unordered_map<ElementType, std::vector<std::string>> type_map = {
        {ElementType::AutomaticConstraint, {"autoconstrain", "autoconstraints"}},
        {ElementType::BodySelect, {"select", "selector"}},
        {ElementType::Constraint, {"constrain", "constraint"}},
        {ElementType::Copy, {"copy", "copy_body"}},
        {ElementType::EveryNStep, {"every"}},
        {ElementType::LoadElement, {"load", "open"}},
        {ElementType::Log, {"log"}},
        {ElementType::LoopBegin, {"loop"}},
        {ElementType::LoopEnd, {"end"}},
        {ElementType::Message, {"print"}},
        {ElementType::OnImprovement, {"on_improvement"}},
        {ElementType::OptimizeStep, {"optimize_step", "optimize_once"}},
        {ElementType::OutputFolder, {"output", "output_folder"}},
        {ElementType::OverlapStrength, {"overlap_strength"}},
        {ElementType::Parameter, {"parameter", "parameter_generator"}},
        {ElementType::RelativeHydration, {"relative_hydration"}},
        {ElementType::Save, {"save", "write"}},
        {ElementType::Seed, {"seed"}},
        {ElementType::SymmetryElement, {"symmetry"}},
        {ElementType::Transform, {"transform", "transformer"}},
    };
    for (const auto& [type, prefixes] : type_map) {
        for (const auto& prefix : prefixes) {
            if (
                line.starts_with(prefix) && (line.size() == prefix.size() || !std::isalpha(line[prefix.size()]))
            ) {return type;}
        }
    }
    throw except::invalid_argument("SequenceParser::get_type: Unknown element type \"" + std::string(line) + "\".");
}

settings::rigidbody::TransformationStrategyChoice get_transform_strategy(std::string_view line) {
    if (line == "rigid_transform" || line == "rigid") {return settings::rigidbody::TransformationStrategyChoice::RigidTransform;}
    if (line == "single_transform" || line == "single") {return settings::rigidbody::TransformationStrategyChoice::SingleTransform;}
    throw except::invalid_argument("SequenceParser::get_transform_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

enum class ConstraintChoice {
    SpecificAtoms,
    Bond,
    BodyCM,
    BodyCMAttractor,
    BodyCMRepeller,
    Unknown
};

ConstraintChoice get_constraint_choice(std::string_view line) {
    if (line == "bond") {return ConstraintChoice::Bond;}
    if (line == "distance") {return ConstraintChoice::SpecificAtoms;}
    if (line == "cm" || line == "center_mass" || line == "center_of_mass") {return ConstraintChoice::BodyCM;}
    if (line == "attract") {return ConstraintChoice::BodyCMAttractor;}
    if (line == "repel") {return ConstraintChoice::BodyCMRepeller;}
    throw except::invalid_argument("SequenceParser::get_constraint_choice: Unknown constraint choice \"" + std::string(line) + "\"");
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Constraint>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    enum class Args {body1, body2, iatom1, iatom2, type, distance};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::body1, {"first", "body1"}},
        {Args::body2, {"second", "body2"}},
        {Args::iatom1, {"iatom1", "iatom_1", "i1"}},
        {Args::iatom2, {"iatom2", "iatom_2", "i2"}},
        {Args::type, {"type", "kind"}},
        {Args::distance, {"distance", "dist", "d"}}
    };
    if (args.size() < 3) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"constraint\". Expected at least 3, but got " + std::to_string(args.size()) + ".");}

    auto body1 = get_arg<std::string>(valid_args[Args::body1], args);
    auto body2 = get_arg<std::string>(valid_args[Args::body2], args);
    auto type  = get_arg<std::string>(valid_args[Args::type], args);
    auto iatom1 = get_arg<int>(valid_args[Args::iatom1], args);
    auto iatom2 = get_arg<int>(valid_args[Args::iatom2], args);
    auto distance = get_arg<double>(valid_args[Args::distance], args);

    ConstraintChoice type_enum = ConstraintChoice::Unknown;
    if (type.found) {type_enum = get_constraint_choice(type.value);}
    {   // input validation
        if (!body1.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"body1\".");}
        if (!body2.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"body2\".");}
        if (iatom1.found && !iatom2.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Argument \"iatom1\" provided without \"iatom2\".");}
        if (!iatom1.found && iatom2.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Argument \"iatom2\" provided without \"iatom1\".");}
        if (type.found && (iatom1.found || iatom2.found)) {throw except::invalid_argument("SequenceParser::parse_arguments: Argument \"type\" cannot be provided together with \"iatom1\" or \"iatom2\".");}
        if (distance.found) {
            if (!type.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Argument \"distance\" cannot be provided without \"type\".");}
        }
        if (type.found) {
            switch (type_enum) {
                case ConstraintChoice::BodyCMAttractor:
                case ConstraintChoice::BodyCMRepeller:
                    if (!distance.found) {
                        throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"distance\" for constraint of type \"" + type.value + "\".");
                    }
                default:
                    break;
            }
        }
    }

    if (!(iatom1.found && iatom2.found)) {
        assert(body1.found && body2.found);
        std::unique_ptr<constraints::Constraint> constraint;
        switch (type_enum) {
            case ConstraintChoice::Bond:
                constraint = factory::create_constraint_bond(
                    &static_cast<Sequencer*>(loop_stack.front())->_get_rigidbody()->molecule,
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body1.value),
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body2.value)
                );
                break;

            case ConstraintChoice::BodyCM:
                constraint = factory::create_constraint_cm(
                    &static_cast<Sequencer*>(loop_stack.front())->_get_rigidbody()->molecule,
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body1.value),
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body2.value)
                );
                break;

            case ConstraintChoice::BodyCMAttractor:
                constraint = factory::create_constraint_attractor(
                    &static_cast<Sequencer*>(loop_stack.front())->_get_rigidbody()->molecule,
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body1.value),
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body2.value),
                    distance.value
                );
                break;

            case ConstraintChoice::BodyCMRepeller:
                constraint = factory::create_constraint_repeller(
                    &static_cast<Sequencer*>(loop_stack.front())->_get_rigidbody()->molecule,
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body1.value),
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body2.value),
                    distance.value
                );
                break;

            case ConstraintChoice::SpecificAtoms:
                constraint = factory::create_constraint(
                    &static_cast<Sequencer*>(loop_stack.front())->_get_rigidbody()->molecule,
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body1.value),
                    loop_stack.back()->_get_sequencer()->setup()._get_body_index(body2.value),
                    iatom1.value,
                    iatom2.value
                );
                break;

            default: 
                throw except::invalid_argument(
                    "SequenceParser::parse_arguments: Unknown constraint type \"" + (type.found ? type.value : "unspecified") + "\". Expected one of {bond, distance, cm, attract, repel}."
                );
        }
        loop_stack.front()->_get_rigidbody()->constraints->add_constraint(std::move(constraint));
    } else {
        assert(body1.found && body2.found && iatom1.found && iatom2.found);
        loop_stack.front()->_get_rigidbody()->constraints->add_constraint(
            factory::create_constraint(
                &static_cast<Sequencer*>(loop_stack.front())->_get_rigidbody()->molecule,
                loop_stack.back()->_get_sequencer()->setup()._get_body_index(body1.value),
                loop_stack.back()->_get_sequencer()->setup()._get_body_index(body2.value),
                iatom1.value,
                iatom2.value
            )
        );
    }
    return nullptr;
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Copy>(const std::unordered_map<std::string, std::vector<std::string>>& args) {}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::OutputFolder>(const std::unordered_map<std::string, std::vector<std::string>>& args) {}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::AutomaticConstraint>(const std::unordered_map<std::string, std::vector<std::string>>& args) {}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Seed>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"seed\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    try {
        int seed = std::stoi(args.begin()->second[0]);
        random::set_seed(seed);
        return nullptr;
    } catch (std::exception&) {
        throw except::invalid_argument("SequenceParser::parse_arguments: \"" + args.begin()->second[0] + "\" cannot be interpreted as an integer seed value.");
    }
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Message>(const std::unordered_map<std::string, std::vector<std::string>>& args) {}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Log>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"log\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    auto message = args.begin()->second[0];
    return std::make_unique<MessageElement>(static_cast<Sequencer*>(loop_stack.front()), message, true);
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::LoadElement>(const std::unordered_map<std::string, std::vector<std::string>>& args) {}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::SymmetryElement>(const std::unordered_map<std::string, std::vector<std::string>>& args) {}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::LoopBegin>(const std::unordered_map<std::string, std::vector<std::string>>& args) {}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Parameter>(const std::unordered_map<std::string, std::vector<std::string>>& args) {}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::BodySelect>(const std::unordered_map<std::string, std::vector<std::string>>& args) {}

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

    // handle anonymous inline: "relative_hydration b1 high" → args["anonymous"] = {"b1", "high"}
    if (args.size() == 1 && args.begin()->first == "anonymous") {
        const auto& vals = args.begin()->second;
        if (vals.size() != 2) {
            throw except::invalid_argument(
                "Element \"relative_hydration\": Inline form expects exactly 2 arguments (body name and level), but got " + std::to_string(vals.size()) + "."
            );
        }
        if (rigidbody->molecule.size_body() != 1) {
            throw except::invalid_argument(
                "Element \"relative_hydration\": Inline form can only be used when there is exactly one body in the molecule. "
                "Use a brace block for multiple bodies."
            );
        }
        return std::make_unique<RelativeHydrationElement>(static_cast<Sequencer*>(loop_stack.front()), std::vector<std::string>{vals[0]}, std::vector<double>{to_value(options[vals[1]])});
    }

    if (args.size() != rigidbody->molecule.size_body()) {
        throw except::invalid_argument(
            "SequenceParser::parse_arguments: Invalid number of arguments for \"relative_hydration\". "
            "Expected " + std::to_string(rigidbody->molecule.size_body()) + ", but got " + std::to_string(args.size()) + "."
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
    observer_ptr<OptimizeStepElement> optimize_step = nullptr;
    if (loop_stack.empty() || !(optimize_step = dynamic_cast<OptimizeStepElement*>(loop_stack.back()))) {
        throw except::invalid_argument("SequenceParser::parse_arguments: \"on_improvement\" must be inside an \"optimize_step\" block.");
    }
    return std::make_unique<OnImprovementElement>(optimize_step);
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Save>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"save\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<SaveElement>(loop_stack.back(), settings::general::output + args.begin()->second[0]);
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::OverlapStrength>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    enum class Args {scaling, distance};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::scaling, {"scaling", "factor"}},
        {Args::distance, {"max", "max_distance", "distance"}}
    };
    if (args.size() != valid_args.size()) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"overlap_strength\". Expected 2, but got " + std::to_string(args.size()) + ".");}

    auto scaling  = get_arg<double>(valid_args[Args::scaling], args);
    auto distance = get_arg<double>(valid_args[Args::distance], args);
    if (!scaling.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"scaling\".");}
    if (!distance.found) {throw except::invalid_argument("SequenceParser::parse_arguments: Missing required argument \"distance\".");}

    static_cast<Sequencer*>(loop_stack.front())->setup().set_overlap_function([a=scaling.value, d=distance.value] (double x) {return x < d ? a*std::pow((d-x)/d, 2) : 0;});
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
    sequencer->setup()._set_config_folder(config.directory());

    std::string line;
    while(!in.eof()) {
        // read the next line
        std::getline(in, line);
        if (line.empty()) {continue;}

        // skip empty lines & comments
        std::unordered_map<std::string, std::vector<std::string>> args;
        
        // tokenize line respecting quoted strings
        static auto tokenize = [] (const std::string& s) -> std::vector<std::string> {
            std::vector<std::string> result;
            std::string current;
            char in_quote = 0;
            bool in_token = false;
            
            for (size_t i = 0; i < s.size(); ++i) {
                char c = s[i];
                
                // handle quotes
                if ((c == '"' || c == '\'') && (i == 0 || s[i-1] != '\\')) {
                    if (in_quote == 0) {
                        in_quote = c;
                        in_token = true;
                    } else if (in_quote == c) {
                        in_quote = 0;
                    } else {
                        current += c;
                    }
                    continue;
                }
                
                // inside quotes: add everything
                if (in_quote != 0) {
                    current += c;
                    continue;
                }
                
                // outside quotes: whitespace delimits tokens
                if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
                    if (in_token) {
                        result.push_back(current);
                        current.clear();
                        in_token = false;
                    }
                    continue;
                }
                
                // regular character
                current += c;
                in_token = true;
            }
            
            if (in_token) {
                result.push_back(current);
            }
            
            return result;
        };
        
        auto tokens = tokenize(line);
        if (tokens.empty() || tokens.front()[0] == '#') {continue;}

        static auto is_opening_brace = [] (char c) {return c == '{' || c == '[' || c == '(';};
        static auto is_closing_brace = [] (char c) {return c == '}' || c == ']' || c == ')';};

        // check if the argument list spans over multiple lines by inspecting the last non-whitespace char
        static auto last_non_ws = [] (const std::string& s) -> char {
            for (int i = static_cast<int>(s.size()) - 1; i >= 0; --i) {
                char c = s[i];
                if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {continue;}
                return c;
            }
            return 0;
        };

        // if the last unquoted non-whitespace character is an opening brace, read a multiline arg block
        if (is_opening_brace(last_non_ws(line))) {
            std::string argline;
            std::getline(in, argline);
            while (!is_closing_brace(last_non_ws(argline))) {
                if (!argline.empty()) {
                    auto sub_tokens = tokenize(argline);
                    if (!sub_tokens.empty()) {
                        int end = static_cast<int>(sub_tokens.size());
                        for (unsigned int i = 0; i < sub_tokens.size(); i++) {
                            if (!sub_tokens[i].empty() && sub_tokens[i][0] == '#') {
                                end = static_cast<int>(i);
                                break;
                            }
                        }
                        args[sub_tokens[0]] = std::vector<std::string>(sub_tokens.begin()+1, sub_tokens.begin()+end);
                    }
                }

                if (in.eof()) { throw except::io_error("SequenceParser::parse: Unescaped argument list starting in line \"" + line + "\"."); }
                std::getline(in, argline);
            }
        }

        // else: all remaining inline tokens are positional anonymous arguments
        // (e.g. "loop 100", "loop copy l1", "copy b2 b1")
        else if (tokens.size() >= 2) {
            args["anonymous"] = std::vector<std::string>(tokens.begin()+1, tokens.end());
        }

        logging::log("Parsed script:");
        std::string parsed_script = tokens[0] + " " + (args.empty() ? "" : "{\n");
        for (const auto& [key, value] : args) {
            parsed_script += "\t" + key + ": ";
            for (const auto& v : value) {
                parsed_script += "\"" + v + "\" ";}
                parsed_script += "\n";
            }
            parsed_script += (args.empty() ? "" : "}");
        logging::log(parsed_script);

        switch (get_type(tokens[0])) {
            case ElementType::Constraint:
                parse_arguments<ElementType::Constraint>(args);
                break;

            case ElementType::AutomaticConstraint:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::AutomaticConstraint>(args));
                break;

            case ElementType::LoopBegin:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::LoopBegin>(args));
                if (auto* loop = dynamic_cast<LoopElement*>(loop_stack.back()->_get_elements().back().get())) {
                    loop_stack.emplace_back(loop); // the returned element _may_ be a CopyLoopElement, which is not modifiable
                }
                break;

            case ElementType::LoopEnd:
                parse_arguments<ElementType::LoopEnd>(args);
                if (loop_stack.size() == 1) {throw except::invalid_argument("Element \"end\": Too many \"end\" statements.");}
                loop_stack.pop_back();
                break;

            case ElementType::Parameter:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::Parameter>(args));
                last_parameter_element = dynamic_cast<ParameterElement*>(loop_stack.back()->_get_elements().back().get());
                break;

            case ElementType::BodySelect:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::BodySelect>(args));
                break;

            case ElementType::Copy:
                parse_arguments<ElementType::Copy>(args);
                break;

            case ElementType::Transform:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::Transform>(args));
                break;

            case ElementType::OptimizeStep:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::OptimizeStep>(args));
                loop_stack.emplace_back(static_cast<LoopElement*>(loop_stack.back()->_get_elements().back().get()));
                break;

            case ElementType::EveryNStep:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::EveryNStep>(args));
                loop_stack.emplace_back(static_cast<LoopElement*>(loop_stack.back()->_get_elements().back().get()));
                break;

            case ElementType::OnImprovement:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::OnImprovement>(args));
                loop_stack.emplace_back(static_cast<LoopElement*>(loop_stack.back()->_get_elements().back().get()));
                break;

            case ElementType::Save:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::Save>(args));
                break;

            case ElementType::LoadElement:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::LoadElement>(args));
                break;

            case ElementType::SymmetryElement:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::SymmetryElement>(args));
                break;

            case ElementType::OverlapStrength:
                parse_arguments<ElementType::OverlapStrength>(args);
                break;

            case ElementType::RelativeHydration:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::RelativeHydration>(args));
                break;

            case ElementType::OutputFolder:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::OutputFolder>(args));
                break;
            
            case ElementType::Message:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::Message>(args));
                break;

            case ElementType::Log:
                loop_stack.back()->_get_elements().emplace_back(parse_arguments<ElementType::Log>(args));
                break;

            case ElementType::Seed:
                parse_arguments<ElementType::Seed>(args);
                break;

            default:
                throw except::invalid_argument("SequenceParser::constructor: Unknown element type \"" + tokens[0] + "\".");
        }
    }
    if (loop_stack.size() != 1) {throw except::invalid_argument("SequenceParser::constructor: Missing \"end\" statements.");}
    return sequencer;
}