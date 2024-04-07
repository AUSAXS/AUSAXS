/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
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

enum class ElementType {
    LoopBegin,
    LoopEnd,
    Parameter,
    BodySelect,
    Transform,
    OptimizeStep,
    EveryNStep,
    Save
};

ElementType get_type(std::string_view line) {
    if (line.starts_with("loop")) {return ElementType::LoopBegin;}
    if (line.starts_with("end")) {return ElementType::LoopEnd;}
    if (line.starts_with("parameter")) {return ElementType::Parameter;}
    if (line.starts_with("body_select")) {return ElementType::BodySelect;}
    if (line.starts_with("transform")) {return ElementType::Transform;}
    if (line.starts_with("optimize_step")) {return ElementType::OptimizeStep;}
    if (line.starts_with("every_n_step")) {return ElementType::EveryNStep;}
    if (line.starts_with("save")) {return ElementType::Save;}
    throw except::invalid_argument("SequenceParser::get_type: Unknown element type \"" + std::string(line) + "\"");
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
    if (line == "both") {return settings::rigidbody::ParameterGenerationStrategyChoice::Simple;}
    throw except::invalid_argument("SequenceParser::get_strategy: Unknown strategy \"" + std::string(line) + "\"");
}

Sequencer SequenceParser::parse(const io::ExistingFile& file) {
    std::ifstream in(file.path());
    if (!in.is_open()) {throw except::io_error("SequenceParser::parse: Could not open file \"" + file.path() + "\".");}

    RigidBody rigidbody = rigidbody::BodySplitter::split("dummy", {});
    Sequencer sequencer("dummy", &rigidbody);

    std::vector<observer_ptr<LoopElement>> loop_stack = {&sequencer};

    std::string line;
    while(std::getline(in, line)) {
        auto tokens = utility::split(line, ' ');
        switch (get_type(tokens[0])) {
            case ElementType::LoopBegin:
                // format: loop <iterations>
                loop_stack.back()->_get_elements().push_back(std::make_unique<LoopElement>(loop_stack.back(), std::stoi(tokens[1])));
                loop_stack.push_back(static_cast<LoopElement*>(loop_stack.back()->_get_elements().back().get()));
                break;

            case ElementType::LoopEnd:
                // format: end
                loop_stack.pop_back();
                break;

            case ElementType::Parameter:
                // format: parameter <iterations> <angstroms> <radians> <strategy> <decay_strategy>
                loop_stack.back()->_get_elements().push_back(
                    std::make_unique<ParameterElement>(
                        loop_stack.back(), 
                        rigidbody::factory::create_parameter_strategy(
                            rigidbody::factory::create_decay_strategy(
                                std::stoi(tokens[1]),
                                get_decay_strategy(tokens[5])
                            ),
                            std::stoi(tokens[2]), 
                            std::stod(tokens[3]),
                            get_parameter_strategy(tokens[4])
                        )
                    )
                );
                break;

            case ElementType::BodySelect:
                // format: body_select <strategy>
                loop_stack.back()->_get_elements().push_back(
                    std::make_unique<BodySelectElement>(
                        loop_stack.back(),
                        rigidbody::factory::create_selection_strategy(
                            sequencer._get_rigidbody(),
                            get_body_select_strategy(tokens[1])
                        )
                    )
                );
                break;

            default:
                throw except::invalid_argument("SequenceParser::constructor: Unknown element type.");
        }
    }
}