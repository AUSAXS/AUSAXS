// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/AdditionalElements.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/All.h>
#include <utility/observer_ptr.h>
#include <utility/Logging.h>
#include <utility/Exceptions.h>
#include <io/ExistingFile.h>

#include <fstream>
#include <unordered_map>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

using ausaxs::except::io_error;
using ausaxs::rigidbody::sequencer::except::parse_error;

namespace {
    // tokenize a line respecting quoted strings; '#' outside quotes starts a comment and discards the rest
    std::vector<std::string> tokenize(const std::string& s) {
        std::vector<std::string> result;
        std::string current;
        char in_quote = 0;
        bool in_token = false;

        for (size_t i = 0; i < s.size(); ++i) {
            char c = s[i];

            // unquoted '#' starts a comment - stop tokenizing
            if (c == '#' && in_quote == 0) {break;}

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
    }

    bool is_opening_brace(char c) {return c == '{' || c == '[' || c == '(';}
    bool is_closing_brace(char c) {return c == '}' || c == ']' || c == ')';}

    // check if the last token is a standalone opening/closing brace (not a value that happens to end with one)
    bool ends_with_opening_brace(const std::vector<std::string>& tokens) {
        if (tokens.empty()) {return false;}
        const auto& t = tokens.back();
        return t.size() == 1 && is_opening_brace(t[0]);
    }

    bool ends_with_closing_brace(const std::vector<std::string>& tokens) {
        if (tokens.empty()) {return false;}
        const auto& t = tokens.back();
        return t.size() == 1 && is_closing_brace(t[0]);
    }
}

enum class ElementType {
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
    throw parse_error("base", "Unknown element \"" + std::string(line) + "\".");
}

std::unique_ptr<Sequencer> SequenceParser::parse(const io::ExistingFile& config) {
    std::ifstream in(config.path());
    if (!in.is_open()) {throw ausaxs::except::io_error("SequenceParser::parse: Could not open file \"" + config.path() + "\".");}

    std::unique_ptr<Sequencer> sequencer = std::make_unique<Sequencer>(); // the main sequencer object

    // the top element of this stack is the current loop element which new elements will be added to
    // note that the sequencer itself is just a dummy loop element with an iteration count of 1
    loop_stack = {sequencer.get()};
    sequencer->setup()._set_config_folder(config.directory());

    std::string line;
    int line_no = 0;
    while(!in.eof()) {
        std::getline(in, line);
        ++line_no;
        if (line.empty()) {continue;}

        ParsedArgs pargs;
        auto tokens = tokenize(line);
        if (tokens.empty()) {continue;}

        if (ends_with_opening_brace(tokens)) { // if the last token ends with an opening brace, read a multiline arg block
            std::string argline;
            std::getline(in, argline);
            ++line_no;
            auto sub_tokens = tokenize(argline);
            while (!ends_with_closing_brace(sub_tokens)) {
                if (!sub_tokens.empty()) {
                    ParsedArgs::Args sub_args;
                    sub_args.line_number = line_no;
                    for (size_t i = 1; i < sub_tokens.size(); ++i) {
                        sub_args.args.push_back({line_no, sub_tokens[i]});
                    }
                    pargs.named[sub_tokens[0]] = std::move(sub_args);
                }

                if (in.eof()) { throw io_error("SequenceParser::parse: Unescaped argument list starting in line \"" + line + "\"."); }
                std::getline(in, argline);
                ++line_no;
                sub_tokens = tokenize(argline);
            }
        }

        // else: all remaining inline tokens are positional anonymous arguments
        else if (tokens.size() >= 2) {
            pargs.inlined.line_number = line_no;
            pargs.inlined.values = std::vector<std::string>(tokens.begin()+1, tokens.end());
        }

        logging::log("Parsed script:");
        std::string parsed_script = tokens[0] + " " + (pargs.named.empty() && pargs.inlined.empty() ? "" : "{\n");
        for (const auto& [key, value] : pargs.named) {
            parsed_script += "\t" + key + ": ";
            for (const auto& v : value.args) {
                parsed_script += "\"" + v.str + "\" ";}
                parsed_script += "\n";
            }
        if (!pargs.inlined.empty()) {
            parsed_script += "\tinline: ";
            for (const auto& v : pargs.inlined.values) {
                parsed_script += "\"" + v + "\" ";}
                parsed_script += "\n";
            }
            parsed_script += (pargs.named.empty() && pargs.inlined.empty() ? "" : "}");
        logging::log(parsed_script);

        auto& elements = loop_stack.back()->_get_elements();
        auto type = get_type(tokens[0]);

        // dispatch maps for the common cases
        using ElementParser = std::unique_ptr<GenericElement>(*)(observer_ptr<LoopElement>, ParsedArgs&&);
        using VoidParser = void(*)(observer_ptr<LoopElement>, ParsedArgs&&);
        static const std::unordered_map<ElementType, ElementParser> element_parsers = {
            {ElementType::AutomaticConstraint, AutoConstraintsElement::_parse},
            {ElementType::BodySelect,          BodySelectElement::_parse},
            {ElementType::Copy,                CopyBodyElement::_parse},
            {ElementType::EveryNStep,           EveryNStepElement::_parse},
            {ElementType::LoadElement,         LoadElement::_parse},
            {ElementType::LoopBegin,           LoopElement::_parse},
            {ElementType::Message,             MessageElement::_parse},
            {ElementType::OnImprovement,       OnImprovementElement::_parse},
            {ElementType::OptimizeStep,        OptimizeStepElement::_parse},
            {ElementType::OutputFolder,        OutputFolderElement::_parse},
            {ElementType::Parameter,           ParameterElement::_parse},
            {ElementType::RelativeHydration,   RelativeHydrationElement::_parse},
            {ElementType::Save,                SaveElement::_parse},
            {ElementType::SymmetryElement,     SymmetryElement::_parse},
            {ElementType::Transform,           TransformElement::_parse},
            {ElementType::Log,                 detail::LogElement::_parse},
        };
        static const std::unordered_map<ElementType, VoidParser> void_parsers = {
            {ElementType::Constraint,          ConstraintElement::_parse},
            {ElementType::LoopEnd,             detail::LoopEndElement::_parse},
            {ElementType::OverlapStrength,     detail::OverlapStrengthElement::_parse},
            {ElementType::Seed,                detail::SeedElement::_parse},
        };

        if (element_parsers.contains(type)) {elements.emplace_back(element_parsers.at(type)(loop_stack.back(), std::move(pargs)));}
        else if (void_parsers.contains(type)) {void_parsers.at(type)(loop_stack.back(), std::move(pargs));}
        switch (type) { // additional work for elements that require it
            case ElementType::LoopBegin:
                // the returned element _may_ be a CopyLoopElement which does not open a new scope
                if (auto* loop = dynamic_cast<LoopElement*>(elements.back().get())) {loop_stack.emplace_back(loop);}
                break;

            case ElementType::LoopEnd:
                if (loop_stack.size() == 1) {throw except::parse_error("end", "Too many \"end\" statements.");}
                loop_stack.pop_back();
                break;

            case ElementType::OptimizeStep:
            case ElementType::EveryNStep:
            case ElementType::OnImprovement:
                loop_stack.emplace_back(static_cast<LoopElement*>(elements.back().get()));
                break;

            default:
                break;
        }
    }
    if (loop_stack.size() != 1) {throw parse_error("base", "Missing \"end\" statement.");}
    return sequencer;
}