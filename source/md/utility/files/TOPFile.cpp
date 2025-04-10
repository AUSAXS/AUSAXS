/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/utility/files/TOPFile.h>
#include <md/utility/Exceptions.h>
#include <utility/Console.h>

#include <algorithm>
#include <fstream>
#include <vector>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::md;

void TOPFile::include(const ITPFile& itp, const std::string& symbol) {
    if (!exists()) {throw except::io_error("TOPFile: \"" + path() + "\" does not exist.");}
    std::ifstream in(path());

    // copy old contents
    std::string line, contents;
    while (std::getline(in, line)) {
        if (line.find(itp.path()) != std::string::npos) {return;}
        contents += line + "\n";}
    in.close();
    
    // insert restraints
    if (symbol.empty()) {
        contents += "\n\n#include \"" + itp.path() + "\"\n";
    } else {
        contents += "\n\n#ifdef " + symbol + "\n";
        contents += "\t#include \"" + itp.path() + "\"\n";
        contents += "#endif\n";
    }

    // write new contents
    std::ofstream out(path());
    out << contents;
    out.close();
}

void TOPFile::include(const ITPFile& itp, const std::string& symbol, const std::string& section) {
    if (!exists()) {throw except::io_error("TOPFile: \"" + path() + "\" does not exist.");}
    std::ifstream in(path());

    // copy old contents
    std::string line;
    std::vector<std::string> contents;
    bool in_block = false;
    bool inserted = false;
    while (std::getline(in, line)) {
        if (line.find(itp.path()) != std::string::npos) {return;}
        if (line.find(section) != std::string::npos) {
            if (inserted) {throw except::io_error("TOPFile: Multiple lines matching \"" + section + "\" found.");}
            in_block = true;
        }

        if (line.empty() && in_block) {
            in_block = false;
            inserted = true;

            if (symbol.empty()) {
                contents.push_back(
                    "\n"
                    "#include \"" + itp.path() + "\"\n"
                );
            } else {
                contents.push_back(
                    "\n"
                    "#ifdef " + symbol + "\n"
                    "\t#include \"" + itp.path() + "\"\n"
                    "#endif\n"
                );
            }
            continue;
        }
        contents.push_back(line);
    }
    if (!inserted) {throw except::io_error("TOPFile: Could not find section \"" + section + "\" in \"" + path() + "\".");}
    in.close();

    // write new contents
    std::string new_contents;
    for (const auto& line : contents) {new_contents += line + "\n";}
    std::ofstream out(path());
    out << new_contents;
}

void TOPFile::include_new_type(const std::vector<ITPFile>& itps, const std::string& symbol) {
    if (!exists()) {throw except::io_error("TOPFile::include: \"" + path() + "\" does not exist.");}
    if (includes.empty()) {includes = discover_includes();}
    if (itps.size() != includes.size()) {
        console::print_critical("TOPFile::include: Number of includes does not match number of supplied itp files.");
        console::print_text_critical("\tPassed include topology files:");
        for (const auto& itp : itps) {
            console::print_text_critical("\t\t" + itp.path());
        }
        console::print_text_critical("\tCurrent include topology files:");
        for (const auto& include : includes) {
            console::print_text_critical("\t\t" + include.path());
        }
        throw except::io_error("TOPFile::include: Number of includes does not match number of supplied itp files.");
    }

    // extract name from filename (scatt_XXX.itp -> XXX)
    std::vector<std::string> name(itps.size());
    {
        std::string type = "";
        for (unsigned int i = 0; i < itps.size(); i++) {
            auto fname = itps[i].filename();
            auto begin = fname.find_first_of('_');
            auto end = fname.find_last_of('.');
            if (begin == std::string::npos || end == std::string::npos) {
                throw except::io_error("TOPFile::include: Could not extract name from filename \"" + itps[i].path() + "\".");
            }
            if (type.empty()) {
                type = fname.substr(0, begin);
                console::print_text_minor("Adding \"" + type + "\" includes to topology file.");
            }
            else if (type != fname.substr(0, begin)) {
                throw except::io_error("TOPFile::include: Type mismatch. All itp files must have the same type. Detected \"" + type + "\" and \"" + fname.substr(0, begin) + "\".");
            }
            name[i] = fname.substr(begin + 1, end - begin - 1);
            // std::cout << "TOPFile: " << itps[i].path << " -> " << name[i] << std::endl;
        }
    }

    // find chain include in topology file and insert itp file after
    std::ifstream in(path());
    std::string line;
    std::vector<std::string> contents;
    std::vector<unsigned int> inserted(name.size(), 0);
    std::vector<bool> already_present(name.size(), 0);
    std::vector<std::string> relative_paths(itps.size());
    std::transform(itps.begin(), itps.end(), relative_paths.begin(), [this] (const ITPFile& itp) {return relative_path(itp.path());});
    while (std::getline(in, line)) {
        contents.push_back(line);
        for (unsigned int i = 0; i < name.size(); i++) {
            if (line.find("topol_" + name[i]) != std::string::npos) {
                // check if itp file was already inserted
                if (line.find(relative_paths[i]) != std::string::npos) {
                    already_present[i] = true;
                    continue;
                }

                if (symbol.empty()) {
                    contents.push_back("#include \"" + relative_paths[i] + "\"");
                } else {
                    contents.push_back("#ifdef " + symbol);
                    contents.push_back("\t#include \"" + relative_paths[i] + "\"");
                    contents.push_back("#endif");
                }
                inserted[i]++;
            }
        }
    }

    // if all itp files were already present, do nothing
    if (std::all_of(already_present.begin(), already_present.end(), [] (bool b) {return b;})) {
        console::print_text_minor("\tAll includes were already present in topology file.");
        return;
    }

    // validate that all itp files were inserted
    auto sum = std::accumulate(inserted.begin(), inserted.end(), 0);
    if (sum < static_cast<int>(name.size())) {
        console::print_critical("TOPFile::include: Could not find include locations for the following itp files:");
        for (const auto& itp : itps) {
            console::print_text_critical('\t' + itp.path());
        }
        throw except::io_error("TOPFile::include: Could not find include locations for all itp files.");
    } else if (sum > static_cast<int>(name.size())) {
        throw except::io_error("TOPFile::include: Multiple include locations found for itp files.");
    }

    // write new contents
    std::string new_contents;
    for (const auto& line : contents) {new_contents += line + "\n";}
    std::ofstream out(path());
    out << new_contents;
    out.close();
}

std::vector<ITPFile> TOPFile::discover_includes() const {
    if (!exists()) {return {};}
    
    // look for itp files in the same directory
    std::vector<ITPFile> includes_dir;
    for (const auto& entry : directory().files()) {
        // only include itp files that contain "topol" in their name
        auto str = entry.filename();
        if (entry.extension() == ".itp" && str.substr(0, 5) == "topol") {
            includes_dir.push_back(ITPFile(entry.path()));
        }
    }

    // check that they are actually included in the topology file
    std::ifstream in(path());
    std::string line;
    std::vector<ITPFile> includes_file;
    while (std::getline(in, line)) {
        if (line.find("#include") != std::string::npos) {
            auto index = line.find("\"");
            if (index == std::string::npos) {throw except::io_error("TOPFile::discover_includes: \"" + path() + "\" contains an include that is not in quotes.");}
            std::string include = line.substr(index + 1);
            index = include.find("\"");
            if (index == std::string::npos) {throw except::io_error("TOPFile::discover_includes: \"" + path() + "\" contains an include quote which is not terminated.");}
            include = include.substr(0, index);

            // check if the include is in the same directory
            for (const auto& itp : includes_dir) {
                File file(include);
                if (itp.filename() == file.filename()) {
                    includes_file.push_back(itp);
                    break;
                }
            }
        }
    }

    if (includes_file.size() != includes_dir.size()) {
        console::print_warning("TOPFile::discover_includes: \"" + path() + "\" does not include all itp files in the same directory.");
        console::indent();
        console::print_text("Includes in file \"" + path() + "\":");
        console::indent();
        for (const auto& itp : includes_file) {
            console::print_text(itp.path());
        }
        console::unindent();
        console::print_text("Include files in directory \"" + directory().path() + "\":");
        console::indent();
        for (const auto& itp : includes_dir) {
            console::print_text(itp.path());
        }
        console::unindent(2);
    }

    return includes_file;
}

void TOPFile::fix_relative_includes() {
    if (!exists()) {throw except::io_error("TOPFile::fix_relative_includes: \"" + path() + "\" does not exist.");}
    fix_relative_includes(*this);
    for (auto& include : includes) {
        fix_relative_includes(include);
    }
}

void TOPFile::fix_relative_includes(const io::File& path) {
    if (!path.exists()) {throw except::io_error("TOPFile::fix_relative_includes: \"" + path.str() + "\" does not exist.");}
    console::print_text_minor("Fixing relative includes in \"" + path.str() + "\".");
    console::indent();
    std::ifstream in(path);

    // copy old contents
    std::string line;
    std::vector<std::string> contents;
    while (std::getline(in, line)) {
        if (line.find("#include") != std::string::npos) {
            auto index = line.find("\"");
            if (index == std::string::npos) {throw except::io_error("TOPFile::fix_relative_includes: \"" + path.str() + "\" contains an include that is not in quotes.");}
            std::string include = line.substr(index + 1);
            index = include.find("\"");
            if (index == std::string::npos) {throw except::io_error("TOPFile::fix_relative_includes: \"" + path.str() + "\" contains an include quote which is not terminated.");}
            include = include.substr(0, index);

            // non-existing files are likely relative to the GROMACS directory and shouldn't be touched
            if (!io::File(include).exists()) {
                contents.push_back(line);
                continue;
            }

            std::string relpath = path.relative_path(include);
            if (!relpath.empty()) {
                console::print_text_minor("Changed relative path from \"" + include + "\" to \"" + relpath + "\"");
                contents.push_back("\t#include \"" + relpath + "\"");
                continue;
            }
        }
        contents.push_back(line);
    }
    in.close();
    console::unindent();

    // write new contents
    std::string new_contents;
    for (const auto& line : contents) {new_contents += line + "\n";}
    std::ofstream out(path);
    out << new_contents;
    out.close();
}

void TOPFile::extract_single_chain() {
    if (!exists()) {throw except::io_error("TOPFile: \"" + path() + "\" does not exist.");}
    std::ifstream in(path());

    std::string newfile = "topol_Protein_chain_A.itp";

    // copy old contents
    std::string line;
    std::vector<std::string> chain_contents;
    std::vector<std::string> top_contents;
    bool in_chain = false;
    while (std::getline(in, line)) {
        if (line.find("[ moleculetype ]") != std::string::npos) {
            in_chain = true;
        } else if (in_chain && line.find("#endif") != std::string::npos) {
            in_chain = false;
            chain_contents.push_back(line);
            top_contents.push_back(
                "; Include chain topologies\n"
                "#include \"" + newfile + "\""
            );
            continue;
        }

        if (in_chain) {
            chain_contents.push_back(line);
        } else {
            top_contents.push_back(line);
        }
    }
    in.close();

    if (chain_contents.empty()) {return;}
    console::print_text("Structure is a single chain. Extracting chain information from topology file to \"" + newfile + "\".");

    // write new contents
    std::string top;
    for (const auto& line : top_contents) {top += line + "\n";}
    std::ofstream out(path());
    out << top;

    std::string itp;
    for (const auto& line : chain_contents) {itp += line + "\n";}
    std::ofstream out2(directory().str() + "/" + newfile);
    out2 << itp;

    // update includes
    includes.push_back(ITPFile(directory().str() + "/" + newfile));
}

std::string TOPFile::copy(const io::Folder& folder) const {
    if (!exists()) {return "";}
    
    // look for itp files in the same directory
    std::vector<ITPFile> includes_dir;
    for (const auto& entry : directory().files()) {
        if (entry.extension() == ".itp") {
            includes_dir.push_back(ITPFile(entry.path()));
        }
    }

    // copy all files to the new directory
    File::copy(folder);
    for (const auto& itp : includes_dir) {
        itp.copy(folder);
    }

    return folder.path() + "/" + filename();
}