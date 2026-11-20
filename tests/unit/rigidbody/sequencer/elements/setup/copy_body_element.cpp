// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <settings/All.h>

#include <support/temp_file.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

TEST_CASE("CopyBodyElement") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 250;

    test::TempFile config(".conf",
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "copy b1 clone\n"
    );
    auto seq = SequenceParser().parse_file(config);
    auto* rb = seq->_get_rigidbody();
    REQUIRE(rb->molecule.size_body() == 2);
    const auto& source = rb->molecule.get_body(0);
    const auto& copy = rb->molecule.get_body(1);

    // the live copy must sit where its stored parameters place it, not on top of its source
    CHECK(copy.get_cm().distance(rb->conformation->absolute_parameters.parameters[1].translation) < 1e-6);
    CHECK(source.get_cm().distance(copy.get_cm()) > 1);
    CHECK(source.get_uid() != copy.get_uid());
}
