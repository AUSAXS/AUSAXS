// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_rigidbody.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>

using namespace ausaxs;

API void rigidbody_config_run(
    const char* path,
    int* status
) {return execute_with_catch([&]() {
    rigidbody::sequencer::SequenceParser().parse(io::ExistingFile(path))->execute();
}, status);}