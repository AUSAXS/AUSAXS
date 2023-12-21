#pragma once

namespace settings {
    namespace molecule {
        extern bool center;               // Decides if the structure will be centered at origo.
        extern bool use_effective_charge; // Decides whether the charge of the displaced water will be included.
        extern bool throw_on_unknown_atom;// Decides whether an exception will be thrown if an unknown atom is encountered.
    }
}