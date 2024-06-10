#include <md/utility/files/TOPFile.h>
#include <md/utility/Exceptions.h>

#include <algorithm>
#include <fstream>
#include <filesystem>
#include <vector>
#include <iostream>
#include <numeric>

using namespace md;

void TOPFile::include(const ITPFile& itp, const std::string& symbol) {
    if (!std::filesystem::exists(path)) {throw except::io_error("TOPFile: \"" + path + "\" does not exist.");}
    std::ifstream in(path);

    // copy old contents
    std::string line, contents;
    while (std::getline(in, line)) {
        if (line.find(itp.path) != std::string::npos) {return;}
        contents += line + "\n";}
    in.close();
    
    // insert restraints
    if (symbol.empty()) {
        contents += "\n\n#include \"" + itp.path + "\"\n";
    } else {
        contents += "\n\n#ifdef " + symbol + "\n";
        contents += "\t#include \"" + itp.path + "\"\n";
        contents += "#endif\n";
    }

    // write new contents
    std::ofstream out(path);
    out << contents;
    out.close();
}

void TOPFile::include(const ITPFile& itp, const std::string& symbol, const std::string& section) {
    if (!std::filesystem::exists(path)) {throw except::io_error("TOPFile: \"" + path + "\" does not exist.");}
    std::ifstream in(path);

    // copy old contents
    std::string line;
    std::vector<std::string> contents;
    bool in_block = false;
    bool inserted = false;
    while (std::getline(in, line)) {
        if (line.find(itp.path) != std::string::npos) {return;}
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
                    "#include \"" + itp.path + "\"\n"
                );
            } else {
                contents.push_back(
                    "\n"
                    "#ifdef " + symbol + "\n"
                    "\t#include \"" + itp.path + "\"\n"
                    "#endif\n"
                );
            }
            continue;
        }
        contents.push_back(line);
    }
    if (!inserted) {throw except::io_error("TOPFile: Could not find section \"" + section + "\" in \"" + path + "\".");}
    in.close();

    // write new contents
    std::string new_contents;
    for (const auto& line : contents) {new_contents += line + "\n";}
    std::ofstream out(path);
    out << new_contents;
}

void TOPFile::include(const std::vector<ITPFile>& itps, const std::string& symbol) {
    if (!std::filesystem::exists(path)) {throw except::io_error("TOPFile: \"" + path + "\" does not exist.");}
    if (includes.empty()) {discover_includes();}
    if (itps.size() != includes.size()) {
        std::cout << "TOPFile: Number of includes does not match number of supplied itp files." << std::endl;
        std::cout << "itps:";
        for (const auto& itp : itps) {std::cout << " " << itp.path;}
        std::cout << std::endl;
        std::cout << "includes:";
        for (const auto& include : includes) {std::cout << " " << include.path;}
        std::cout << std::endl;
        throw except::io_error("TOPFile: Number of includes does not match number of supplied itp files.");
    }

    // extract name from filename (scatt_XXX.itp -> XXX)
    std::vector<std::string> name(itps.size());
    for (unsigned int i = 0; i < itps.size(); i++) {
        auto fname = itps[i].filename();
        auto begin = fname.find_first_of('_');
        auto end = fname.find_last_of('.');
        if (begin == std::string::npos || end == std::string::npos) {
            throw except::io_error("TOPFile: Could not extract name from filename \"" + itps[i].path + "\".");
        }
        name[i] = fname.substr(begin + 1, end - begin - 1);
        // std::cout << "TOPFile: " << itps[i].path << " -> " << name[i] << std::endl;
    }

    // find chain include in topology file and insert itp file after
    std::ifstream in(path);
    std::string line;
    std::vector<std::string> contents;
    std::vector<unsigned int> inserted(name.size(), 0);
    std::vector<bool> already_present(name.size(), 0);
    std::vector<std::string> relative_paths(itps.size());
    std::transform(itps.begin(), itps.end(), relative_paths.begin(), [this] (const ITPFile& itp) {return relative_path(itp.path);});
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
    if (std::all_of(already_present.begin(), already_present.end(), [] (bool b) {return b;})) {return;}

    // validate that all itp files were inserted
    auto sum = std::accumulate(inserted.begin(), inserted.end(), 0);
    if (sum < static_cast<int>(name.size())) {
        std::cout << "TOPFile: Could not find include locations for the following itp files:" << std::endl;
        for (const auto& itp : itps) {std::cout << itp.path << std::endl;}
        throw except::io_error("TOPFile::include: Could not find include locations for all itp files.");
    } else if (sum > static_cast<int>(name.size())) {
        throw except::io_error("TOPFile::include: Multiple include locations found for itp files.");
    }

    // write new contents
    std::string new_contents;
    for (const auto& line : contents) {new_contents += line + "\n";}
    std::ofstream out(path);
    out << new_contents;
}

std::vector<ITPFile> TOPFile::discover_includes() const {
    if (!exists()) {return {};}
    
    // look for itp files in the same directory
    std::vector<ITPFile> includes_dir;
    for (const auto& entry : std::filesystem::directory_iterator(parent_path())) {
        // only include itp files that contain "topol" in their name
        auto str = entry.path().filename().string();
        if (entry.path().extension() == ".itp" && str.substr(0, 5) == "topol") {
            includes_dir.push_back(ITPFile(entry.path()));
        }
    }

    // check that they are actually included in the topology file
    std::ifstream in(path);
    std::string line;
    std::vector<ITPFile> includes_file;
    while (std::getline(in, line)) {
        if (line.find("#include") != std::string::npos) {
            auto index = line.find("\"");
            if (index == std::string::npos) {throw except::io_error("TOPFile: \"" + path + "\" contains an include that is not in quotes.");}
            std::string include = line.substr(index + 1);
            index = include.find("\"");
            if (index == std::string::npos) {throw except::io_error("TOPFile: \"" + path + "\" contains an include quote which is not terminated.");}
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
        std::cout << "TOPFile: \"" << path << "\" does not include all itp files in the same directory." << std::endl;
        std::cout << "Includes in file \"" << path << "\":" << std::endl;
        for (const auto& itp : includes_file) {
            std::cout << "\t" << itp.path << std::endl;
        }
        std::cout << "Includes in directory \"" << parent_path() << "\":" << std::endl;
        for (const auto& itp : includes_dir) {
            std::cout << "\t" << itp.path << std::endl;
        }
        throw except::io_error("TOPFile: \"" + path + "\" does not include all itp files in the same directory.");
    }

    return includes_file;
}

void TOPFile::fix_relative_includes() {
    fix_relative_includes(*this);
    for (auto& include : includes) {
        fix_relative_includes(include);
    }
}

void TOPFile::fix_relative_includes(const std::string& path) {
    if (!std::filesystem::exists(path)) {throw except::io_error("TOPFile::fix_relative_includes: \"" + path + "\" does not exist.");}
    std::ifstream in(path);

    // copy old contents
    File f(path);
    std::string line;
    std::vector<std::string> contents;
    while (std::getline(in, line)) {
        if (line.find("#include") != std::string::npos) {
            auto index = line.find("\"");
            if (index == std::string::npos) {throw except::io_error("TOPFile::fix_relative_includes: \"" + path + "\" contains an include that is not in quotes.");}
            std::string include = line.substr(index + 1);
            index = include.find("\"");
            if (index == std::string::npos) {throw except::io_error("TOPFile::fix_relative_includes: \"" + path + "\" contains an include quote which is not terminated.");}
            include = include.substr(0, index);

            std::string relpath = f.relative_path(include);
            // std::cout << "Relative path from \"" << path << "\" to \"" << include << "\": \"" << relpath << "\"" << std::endl;
            if (!relpath.empty()) {
                contents.push_back("\t#include \"" + f.relative_path(include) + "\"");
            } else {
                contents.push_back(line);
            }
            continue;
        }
        contents.push_back(line);
    }
    in.close();

    // write new contents
    std::string new_contents;
    for (const auto& line : contents) {new_contents += line + "\n";}
    std::ofstream out(path);
    out << new_contents;
}

void TOPFile::extract_single_chain() {
    if (!exists()) {throw except::io_error("TOPFile: \"" + path + "\" does not exist.");}
    std::ifstream in(path);

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

    // write new contents
    std::string top;
    for (const auto& line : top_contents) {top += line + "\n";}
    std::ofstream out(path);
    out << top;

    std::string itp;
    for (const auto& line : chain_contents) {itp += line + "\n";}
    std::ofstream out2(parent_path() + "/" + newfile);
    out2 << itp;

    // update includes
    includes.push_back(ITPFile(parent_path() + "/" + newfile));
}

std::string TOPFile::copy(const Folder& folder) const {
    if (!exists()) {return "";}
    
    // look for itp files in the same directory
    std::vector<ITPFile> includes_dir;
    for (const auto& entry : std::filesystem::directory_iterator(parent_path())) {
        if (entry.path().extension() == ".itp") {
            includes_dir.push_back(ITPFile(entry.path()));
        }
    }

    // copy all files to the new directory
    File::copy(folder);
    for (const auto& itp : includes_dir) {
        itp.copy(folder);
    }

    return folder.path + "/" + filename();
}