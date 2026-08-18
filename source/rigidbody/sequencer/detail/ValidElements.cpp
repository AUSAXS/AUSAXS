// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/ValidElements.h>
#include <rigidbody/sequencer/detail/AdditionalElements.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/All.h>

#include <algorithm>
#include <map>
#include <cassert>

using namespace ausaxs::rigidbody::sequencer::detail;

namespace {
    std::string quote_join(std::vector<std::string> values) {
        std::sort(values.begin(), values.end()); // the keys arrive from an unordered_map, so sort for a reproducible message
        std::string joined;
        for (const auto& v : values) {
            if (!joined.empty()) {joined += ", ";}
            joined += "\"" + v + "\"";
        }
        return joined;
    }
}

const std::map<ElementType, std::vector<std::string>>& get_type_map() {
    static std::map<ElementType, std::vector<std::string>> type_map = {
        {ElementType::AutomaticConstraint, {"autoconstrain", "autoconstraints"}},
        {ElementType::BodySelect, {"select", "selector"}},
        {ElementType::Constraint, {"constrain", "constraint"}},
        {ElementType::ConvertToSymmetry, {"convert_to_symmetry"}},
        {ElementType::Copy, {"copy", "copy_body"}},
        {ElementType::Delete, {"delete"}},
        {ElementType::EveryNStep, {"every"}},
        {ElementType::LoadElement, {"load", "open"}},
        {ElementType::Log, {"log"}},
        {ElementType::LoopBegin, {"loop"}},
        {ElementType::LoopEnd, {"end"}},
        {ElementType::Merge, {"merge"}},
        {ElementType::Message, {"print"}},
        {ElementType::OnImprovement, {"on_improvement"}},
        {ElementType::OptimizeStep, {"optimize_step", "optimize_once"}},
        {ElementType::OutputFolder, {"output", "output_folder"}},
        {ElementType::OverlapStrength, {"overlap_strength"}},
        {ElementType::Parameter, {"parameter", "parameter_generator"}},
        {ElementType::RelativeHydration, {"relative_hydration"}},
        {ElementType::Rename, {"rename"}},
        {ElementType::Save, {"save", "write"}},
        {ElementType::Seed, {"seed"}},
        {ElementType::Split, {"split"}},
        {ElementType::SymmetryElement, {"symmetry"}},
        {ElementType::Transform, {"transform", "transformer"}},
        {ElementType::Update, {"update"}},
    };
    return type_map;
};

std::vector<std::string> ausaxs::rigidbody::sequencer::detail::valid_elements() {
    const auto& type_map = get_type_map();
    std::vector<std::string> elements; 
    elements.reserve(type_map.size() * 2); // reserve space for all prefixes (most elements have 2)
    for (int i = 0; i < static_cast<int>(ElementType::COUNT); ++i) {
        ElementType type = static_cast<ElementType>(i);
        assert(type_map.contains(type));
        for (const auto& prefix : type_map.at(type)) {
            elements.emplace_back(prefix);
        }
    }
    return elements;
}

std::vector<std::string> ausaxs::rigidbody::sequencer::detail::valid_arguments(ElementType type) {
    switch (type) {
        case ElementType::AutomaticConstraint: return AutoConstraintsElement::_valid_arguments();
        case ElementType::BodySelect:          return BodySelectElement::_valid_arguments();
        case ElementType::Constraint:          return ConstraintElement::_valid_arguments();
        case ElementType::ConvertToSymmetry:   return ConvertToSymmetryElement::_valid_arguments();
        case ElementType::Copy:                return CopyBodyElement::_valid_arguments();
        case ElementType::Delete:              return DeleteElement::_valid_arguments();
        case ElementType::EveryNStep:          return EveryNStepElement::_valid_arguments();
        case ElementType::LoadElement:         return LoadElement::_valid_arguments();
        case ElementType::Log:                 return detail::LogElement::_valid_arguments();
        case ElementType::LoopBegin:           return LoopElement::_valid_arguments();
        case ElementType::LoopEnd:             return detail::LoopEndElement::_valid_arguments();
        case ElementType::Merge:               return MergeElement::_valid_arguments();
        case ElementType::Message:             return MessageElement::_valid_arguments();
        case ElementType::OnImprovement:       return OnImprovementElement::_valid_arguments();
        case ElementType::OptimizeStep:        return OptimizeStepElement::_valid_arguments();
        case ElementType::OutputFolder:        return OutputFolderElement::_valid_arguments();
        case ElementType::OverlapStrength:     return detail::OverlapStrengthElement::_valid_arguments();
        case ElementType::Parameter:           return ParameterElement::_valid_arguments();
        case ElementType::RelativeHydration:   return RelativeHydrationElement::_valid_arguments();
        case ElementType::Rename:              return RenameElement::_valid_arguments();
        case ElementType::Save:                return SaveElement::_valid_arguments();
        case ElementType::Seed:                return detail::SeedElement::_valid_arguments();
        case ElementType::Split:               return SplitElement::_valid_arguments();
        case ElementType::SymmetryElement:     return SymmetryElement::_valid_arguments();
        case ElementType::Transform:           return TransformElement::_valid_arguments();
        case ElementType::Update:              return UpdateElement::_valid_arguments();
        default: return {};
    }
}

void ausaxs::rigidbody::sequencer::detail::validate_named_arguments(ElementType type, std::string_view element, const ParsedArgs& args) {
    if (args.named.empty()) {return;}

    auto valid = valid_arguments(type);
    std::vector<std::string> unknown;
    int line_number = -1;
    for (const auto& [key, value] : args.named) {
        if (std::find(valid.begin(), valid.end(), key) != valid.end()) {continue;}
        unknown.push_back(key);
        if (line_number == -1 || value.line_number < line_number) {line_number = value.line_number;}
    }
    if (unknown.empty()) {return;}

    std::string msg = "Unknown argument" + std::string(1 < unknown.size() ? "s " : " ") + quote_join(std::move(unknown));
    if (0 <= line_number) {msg += " (first at line " + std::to_string(line_number) + ")";}
    msg += valid.empty()
        ? ". This element takes no named arguments."
        : ". Valid arguments are: " + quote_join(std::move(valid)) + ".";
    throw except::parse_error(element, msg);
}

ElementType ausaxs::rigidbody::sequencer::detail::get_type(std::string_view line) {
    const auto& type_map = get_type_map();
    for (const auto& [type, prefixes] : type_map) {
        for (const auto& prefix : prefixes) {
            if (line.starts_with(prefix) && (line.size() == prefix.size() || !std::isalpha(line[prefix.size()]))) {
                return type;
            }
        }
    }
    throw except::parse_error("base", "Unknown element \"" + std::string(line) + "\".");
}