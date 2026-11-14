// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/detail/dataset/DATWriter.h>

#include <dataset/Dataset.h>
#include <io/File.h>
#include <io/detail/dataset/QUnitAnnotation.h>

#include <iomanip>
#include <sstream>

using namespace ausaxs;

void io::detail::dataset::DATWriter::write(const Dataset& data, const io::File& path, const std::string& header) {
    std::stringstream output;

    // annotate with the unit our q-values are in, so a re-read cannot rescale them
    output << qunit::unit_line(header);

    // write header
    if (!header.empty()) {
        output << header << "\n";
    }
    output << "\n";

    // write data
    int N = data.size_rows(), M = data.size_cols();
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M-1; j++) {
            output << std::left << std::setw(16) << std::setprecision(8) << std::scientific << data.index(i, j) << "\t";
        }
        output << data.index(i, M-1) << "\n";
    }
    path.create(output.str());
}
