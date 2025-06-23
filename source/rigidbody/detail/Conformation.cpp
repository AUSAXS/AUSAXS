#include <rigidbody/detail/Conformation.h>

using namespace ausaxs::rigidbody::detail;

Conformation::Conformation() = default;

Conformation::Conformation(observer_ptr<const Rigidbody> rigidbody) 
    : configuration(rigidbody), original_conformation(rigidbody->molecule.get_bodies()) 
{}

Conformation::~Conformation() = default;