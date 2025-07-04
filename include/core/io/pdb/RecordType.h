// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::io::pdb {
    enum class RecordType {HEADER, ATOM, WATER, TERMINATE, FOOTER, NOTYPE};
}