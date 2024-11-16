/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/programs/gmx.h>

#include <fstream>
#include <chrono>
#include <ctime>

using namespace ausaxs;
using namespace ausaxs::md;

shell::Command& gmx::command() {
    validate();
    for (auto& opt : options) {
        cmd.append(opt->get());
    }
    cmd.append("-nocopyright -quiet");
    return cmd;
}

// to encasuplate the entire command with '', we need to escape all ' characters in the command itself
auto escape_single_quotes = [] (std::string_view cmd) {
    std::string escaped_str;
    for (auto c : cmd) {
        if (c == '\'') {
            escaped_str += "'\\''";
            continue;
        }
        escaped_str += c;
    }
    return escaped_str;
};

std::string gmx::execute() {
    auto cmd = command();
    if (log) {
        write_cmdlog(cmd.get());

        // ensure we're running as bash to use 'set -o pipefail' to avoid the piping hiding errors
        // this will crash if using a shell script, so run everything in bash and encapsulate entire cmd with ''
        cmd.set(escape_single_quotes(cmd.get()));
        cmd.prepend("exec bash -c 'set -o pipefail; ");
        cmd.append("2>&1 | tee -a " + outputlog + "'");
    }

    auto result = cmd.execute();
    if (result.exit_code != 0) {
        throw except::io_error("gmx::gmx: Error executing command: \"" + cmd.get() + "\".");
    }
    return result.out;
}

bool gmx::valid_executable() {
    auto tmp = cmd.append("-version");
    write_cmdlog(tmp.get());
    auto res = tmp.execute();
    return res.exit_code == 0;
};

void gmx::set_logfile(const io::File& log, const io::File& cmdlog) {
    if (log.extension() != ".log" || cmdlog.extension() != ".log"){
        throw except::invalid_format("gmx::gmx: Output log file must have extension \".log\".");
    }
    log.remove();
    gmx::outputlog = log;
    gmx::cmdlog = cmdlog;

    log.create();
    cmdlog.create();
    auto time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    auto msg = 
        "\n#################################################"
        "\n   Program started on " + std::string(std::ctime(&time)) + 
        "#################################################";
    write_log(msg);
    write_cmdlog(msg);
    gmx::log = true;

    if (!gmx::outputlog.exists() || !gmx::cmdlog.exists()) {
        throw except::io_error("gmx::gmx: Could not create log files.");
    }
}

void gmx::validate() const {}

void gmx::write_cmdlog(std::string_view entry) {
    if (cmdlog.path().empty()) {return;}
    std::ofstream log(cmdlog, std::ios_base::app);
    log << entry << std::endl;
}

void gmx::write_log(std::string_view entry) {
    if (outputlog.path().empty()) {return;}
    std::ofstream log(outputlog, std::ios_base::app);
    log << entry << std::endl;
}

std::string option::to_string(Forcefield opt) {
    switch (opt) {
        case Forcefield::AMBER99SB:
            return "amber99sb";
        case Forcefield::AMBER99SB_ILDN:
            return "amber99sb-ildn";
        default:
            throw except::unknown_type("gmx::to_string: Unknown forcefield. (Did you forget to add it to the enum?)");
    }
}

std::string option::to_string(WaterModel opt) {
    switch (opt) {
        case WaterModel::TIP3P:
            return "tip3p";
        case WaterModel::TIP4P:
            return "tip4p";
        case WaterModel::TIP4P2005:
            return "tip4p2005";
        default:
            throw except::unknown_type("gmx::to_string: Unknown water model. (Did you forget to add it to the enum?)");
    }
}

std::string option::to_string(BoxType opt) {
    switch (opt) {
        case BoxType::CUBIC:
            return "cubic";
        case BoxType::TRICLINIC:
            return "triclinic";
        case BoxType::DODECAHEDRON:
            return "dodecahedron";
        case BoxType::OCTAHEDRON:
            return "octahedron";
        default:
            throw except::unknown_type("gmx::to_string: Unknown box type. (Did you forget to add it to the enum?)");
    }
}

std::string option::to_string(Cation opt) {
    switch (opt) {
        case Cation::NA:
            return "NA";
        default:
            throw except::unknown_type("gmx::to_string: Unknown cation type. (Did you forget to add it to the enum?)");
    }
}

std::string option::to_string(Anion opt) {
    switch (opt) {
        case Anion::CL:
            return "CL";
        default:
            throw except::unknown_type("gmx::to_string: Unknown anion type. (Did you forget to add it to the enum?)");
    }
}