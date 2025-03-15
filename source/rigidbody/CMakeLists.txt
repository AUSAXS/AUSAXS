add_library(ausaxs_rigidbody OBJECT 
	"RigidBody.cpp"
	"BodySplitter.cpp"

	"constraints/Constraint.cpp"
	"constraints/ConstraintManager.cpp"
	"constraints/DistanceConstraint.cpp"
	"constraints/FixedConstraint.cpp"
	"constraints/LinkConstraint.cpp"
	"constraints/OverlapConstraint.cpp"
	"constraints/generation/ConstraintGenerationFactory.cpp"
	"constraints/generation/ConstraintGenerationStrategy.cpp"
	"constraints/generation/LinearConstraints.cpp"
	"constraints/generation/NoConstraints.cpp"
	"constraints/generation/VolumetricConstraints.cpp"

	"detail/BestConf.cpp"
	
	"parameters/ParameterGenerationFactory.cpp"
	"parameters/ParameterGenerationStrategy.cpp"
	"parameters/decay/DecayFactory.cpp"
	"parameters/decay/ExponentialDecay.cpp"
	"parameters/decay/LinearDecay.cpp"

	"selection/BodySelectFactory.cpp"
	"selection/BodySelectStrategy.cpp"
	"selection/ManualSelect.cpp"
	"selection/RandomBodySelect.cpp"
	"selection/RandomConstraintSelect.cpp"
	"selection/SequentialBodySelect.cpp"
	"selection/SequentialConstraintSelect.cpp"

	"sequencer/BodySelectElement.cpp"
	"sequencer/ConstraintIteratorElement.cpp"
	"sequencer/ConstraintIteratorElementCallback.cpp"
	"sequencer/EveryNStepElement.cpp"
	"sequencer/LoopElement.cpp"
	"sequencer/LoopElementCallback.cpp"
	"sequencer/OnImprovementElement.cpp"
	"sequencer/OptimizeStepElement.cpp"
	"sequencer/ParameterElement.cpp"
	"sequencer/ParameterElementCallback.cpp"
	"sequencer/SaveElement.cpp"
	"sequencer/Sequencer.cpp"
	"sequencer/TransformElement.cpp"
	"sequencer/detail/SequenceParser.cpp"
	"sequencer/setup/AutoConstraintsElement.cpp"
	"sequencer/setup/ConstraintElement.cpp"
	"sequencer/setup/LoadElement.cpp"
	"sequencer/setup/LoadExistingElement.cpp"
	"sequencer/setup/RelativeHydrationElement.cpp"
	"sequencer/setup/SetupElement.cpp"
	"sequencer/setup/SymmetryElement.cpp"
	
	"transform/RigidTransform.cpp"
	"transform/SingleTransform.cpp"
	"transform/TransformFactory.cpp"
	"transform/TransformGroup.cpp"
	"transform/TransformStrategy.cpp"
)

target_include_directories(ausaxs_rigidbody PUBLIC "${CMAKE_SOURCE_DIR}/include/rigidbody")
target_link_libraries(ausaxs_rigidbody ausaxs_core)
set_target_properties(ausaxs_rigidbody PROPERTIES POSITION_INDEPENDENT_CODE 1)
