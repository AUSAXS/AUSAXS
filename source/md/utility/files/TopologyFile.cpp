#include <md/utility/files/TopologyFile.h>
#include <md/utility/Exceptions.h>
#include <md/utility/files/File.h>

#include <filesystem>
#include <fstream>

using namespace gmx;

void TopologyFile::discover_includes() {
    if (!std::filesystem::exists(path)) {return;}
    std::ifstream in(path);

    unsigned int lineno = 0;
    std::string line;
    std::vector<Include> includes;
    while (std::getline(in, line)) {
        if (line.find("#include") != std::string::npos) {
            std::string include = line.substr(line.find("\"") + 1);
            include = include.substr(0, include.find("\""));
            includes.push_back({lineno, include});
        }
        lineno++;
    }
    this->includes = includes;
}

void TopologyFile::fix_relative_includes() {
    if (!std::filesystem::exists(path)) {throw except::io_error("gmx::TopologyFile: \"" + path + "\" does not exist.");}
    if (includes.empty()) {discover_includes();}

    std::ifstream in(path);

    // copy old contents
    std::vector<std::string> contents;
    std::string line;
    while (std::getline(in, line)) {contents.push_back(line);}
    in.close();

    // fix includes
    detail::File file(path);
    for (const auto& include : includes) {
        std::string rpath = file.relative_path(include.path);
        if (!rpath.empty()) {
            contents[include.line] = "\t#include \"" + rpath + "\"";
        }

        // recursively fix includes in included files
        if (include.path.find("topol") != std::string::npos) {
            TopologyFile include_file(include.path);
            include_file.fix_relative_includes();
        }
    }

    // write new contents
    std::string new_contents;
    for (const auto& line : contents) {new_contents += line + "\n";}
    std::ofstream out(path);
    out << new_contents;
    out.close();
}

// void TopologyFile::replace_restraints(const std::vector<std::string>& restraints) {
//     if (!std::filesystem::exists(path)) {throw except::io_error("gmx::TopologyFile: \"" + path + "\" does not exist.");}
//     std::ifstream in(path);

//     // copy old contents
//     std::vector<std::string> contents;
//     std::string line;
//     while (std::getline(in, line)) {contents.push_back(line);}
//     in.close();

//     // recursively replace restraints in included files
//     if (restraints.size() > 1) {
//         unsigned int count = 0;
//         for (const auto& include : includes) {
//             if (include.path.find("topol") != std::string::npos) {count++;}
//         }
//         if (count != restraints.size()) {throw except::io_error("gmx::TopologyFile: number of restraints does not match number of topology files.");}        

//         count = 0;
//         for (const auto& include : includes) {
//             if (include.path.find("topol") != std::string::npos) {
//                 TopologyFile include_file(include.path);
//                 include_file.replace_restraints({restraints[count++]});
//             }
//         }
//         return;
//     }

//     // we only have a single restraint file

//     // check if restraints are already included
//     for (const auto& include : includes) {
//         if (include.path == restraints) {return;}
//     }

//     // replace old restraints with new restraints
//     for (const auto& include : includes) {
//         if (include.path.find("posre") != std::string::npos) {
//             if (line.size() > 15) {continue;}
//             contents[include.line-1] = "#ifdef POSRESBACKBONE";
//             contents[include.line]   = "\t#include \"" + restraints + "\"";
//             contents[include.line+1] = "#endif";
//         }
//     }

//     // copy old contents
//     std::string line;
//     std::vector<std::string> contents;
//     unsigned int posres = 0;
//     while (std::getline(in, line)) {
//         if (line.find("#ifdef POSRESBACKBONE") != std::string::npos) {
//             return;
//         }
//         if (line.find("#ifdef POSRES") != std::string::npos) {
//             if (line.size() < 16) { // ignore #ifdef POSRES_*
//                 posres = contents.size();
//             }
//         }
//         contents.push_back(line);
//     }
//     in.close();
    
//     // determine where to insert restraints
//     if (posres != 0) {
//         contents[posres] =   "#ifdef POSRESBACKBONE";
//         contents[posres+1] = "\t#include \"" + restraints.path + "\"";
//         contents[posres+2] = "#endif";
//     } else {
//         throw except::io_error("gmx::TOPFile::insert_restraints: Failed.");
//     }

//     // write new contents
//     std::string new_contents;
//     for (const auto& line : contents) {new_contents += line + "\n";}
//     std::ofstream out(path);
//     out << new_contents;
//     out.close();
// }
