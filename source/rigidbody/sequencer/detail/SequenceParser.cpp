// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
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

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::EveryNStep>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"every_n_step\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<EveryNStepElement>(loop_stack.back(), std::stoi(args.begin()->second[0]));
}

template<>
std::unique_ptr<GenericElement> SequenceParser::parse_arguments<ElementType::Save>(const std::unordered_map<std::string, std::vector<std::string>>& args) {
    if (args.size() != 1) {throw except::invalid_argument("SequenceParser::parse_arguments: Invalid number of arguments for \"save\". Expected 1, but got " + std::to_string(args.size()) + ".");}
    return std::make_unique<SaveElement>(loop_stack.back(), settings::general::output + args.begin()->second[0]);
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