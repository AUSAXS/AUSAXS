// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// The benchmarks of the AUSAXS paper: run time of a full SAS fit against CRYSOL, FoXS and Pepsi-SAXS, and run time and
// accuracy of the pure Debye sum against DebyeCalculator and AES-Debye. See tests/benchmark/paper/README.md for the setup
// of the third-party programs and the input bundle.
//
// Run from the repository root:
//     cmake --build build-release --target benchmark_paper -j8
//     ./build-release/tests/benchmark/bin/benchmark_paper "[paper]"
//     python3 tests/benchmark/plot_paper.py output/paper_benchmark
//
// Unlike benchmark_pr, this does not use Catch2's BENCHMARK macro: most of what it times are other programs, run as
// whole processes or in their own Python interpreter, and all of them must be reported in one format. Every
// measurement is one warm-up run followed by AUSAXS_BENCH_REPS timed runs, all of which are written to results.jsonl,
// so the plot script can choose the statistic. The file is appended to and read back on start, so an interrupted run
// resumes where it stopped. Configuration is through environment variables, listed in Config below.

#include <catch2/catch_test_macros.hpp>

#include <data/Molecule.h>
#include <dataset/SimpleDataset.h>
#include <fitter/FitResult.h>
#include <fitter/SmartFitter.h>
#include <gpu/GPULoader.h>
#include <hist/scattering_model/DebyeGradient.h>
#include <hist/scattering_model/IScatteringModel.h>
#include <io/ExistingFile.h>
#include <io/File.h>
#include <settings/All.h>

#include <spawn.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>  // environ

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>


using namespace ausaxs;
namespace fs = std::filesystem;

namespace {
    //###################//
    //### CONFIGURATION ###//
    //###################//

    std::string env(const char* name, std::string fallback = "") {
        const char* v = std::getenv(name);
        return v && *v ? std::string(v) : fallback;
    }

    std::vector<double> env_list(const char* name, const std::vector<double>& fallback) {
        auto s = env(name);
        if (s.empty()) {return fallback;}
        std::vector<double> out;
        std::stringstream ss(s);
        for (std::string item; std::getline(ss, item, ',');) {out.push_back(std::stod(item));}
        return out;
    }

    // An external program: its executable, and an optional directory to prepend to its LD_LIBRARY_PATH, for old
    // binaries whose shared libraries the system no longer ships.
    struct Tool {
        std::string exe;
        std::string libs;
        explicit operator bool() const {return !exe.empty();}
    };

    struct Config {
        fs::path manifest   = env("AUSAXS_BENCH_MANIFEST", "tests/benchmark/paper/manifest.tsv"); // the SAS structures
        fs::path data       = env("AUSAXS_BENCH_DATA", "data");                       // the root the manifest's paths are relative to
        fs::path out        = env("AUSAXS_BENCH_OUT", "output/paper_benchmark");         // results.jsonl, meta.json, iq/, work/
        std::string ausaxs  = env("AUSAXS_BENCH_AUSAXS", AUSAXS_BENCH_DEFAULT_CLI);      // the CLI, for the end-to-end series
        std::string python  = env("AUSAXS_BENCH_PYTHON", "python3");                    // with debyecalculator and/or aesdebye
        std::string adapter = env("AUSAXS_BENCH_ADAPTER", "tests/benchmark/paper/competitor.py");
        int reps            = std::stoi(env("AUSAXS_BENCH_REPS", "5"));
        int threads         = std::stoi(env("AUSAXS_BENCH_THREADS", std::to_string(std::thread::hardware_concurrency())));
        double budget       = std::stod(env("AUSAXS_BENCH_BUDGET", "300"));  // s; a single run past this ends the series
        bool gpu            = env("AUSAXS_BENCH_GPU", "1") != "0";            // also run the GPU variants
        std::string device  = env("AUSAXS_BENCH_DEVICE", "cuda");             // torch/AES device for the competitors' GPU runs

        Tool pepsi  {env("AUSAXS_BENCH_PEPSI"),  env("AUSAXS_BENCH_PEPSI_LIBS")};
        Tool foxs   {env("AUSAXS_BENCH_FOXS"),   env("AUSAXS_BENCH_FOXS_LIBS")};
        Tool crysol {env("AUSAXS_BENCH_CRYSOL"), env("AUSAXS_BENCH_CRYSOL_LIBS")};
        std::set<std::string> competitors = [] {  // which Python libraries to run: dc, aes
            std::set<std::string> s;
            std::stringstream ss(env("AUSAXS_BENCH_COMPETITORS", "dc,aes"));
            for (std::string item; std::getline(ss, item, ',');) {if (!item.empty()) {s.insert(item);}}
            return s;
        }();

        // Radii in Å of the fcc silver spheres of the pure Debye series. 10 Å is ~250 atoms; 93 Å is ~200k.
        std::vector<double> radii          = env_list("AUSAXS_BENCH_RADII", {10, 15, 20, 30, 40, 50, 58, 75, 93, 120, 150});
        std::vector<double> accuracy_radii = env_list("AUSAXS_BENCH_ACCURACY_RADII", {10, 20, 30});
        std::vector<double> bin_widths     = env_list("AUSAXS_BENCH_BIN_WIDTHS", {0.5, 0.25, 0.1, 0.05, 0.02, 0.01});
    };

    const Config& config() {
        static Config c;
        return c;
    }

    // The q grid of the pure Debye series, shared by all three programs: q = 0.01 i for i = 1..1000, i.e. up to 10 Å⁻¹.
    // AES-Debye builds it as start + i (end - start)/steps with the end excluded, DebyeCalculator with arange.
    constexpr int n_q = 1000;
    constexpr double q_step = 0.01;
    std::vector<double> debye_q() {
        std::vector<double> q(n_q);
        for (int i = 0; i < n_q; ++i) {q[i] = q_step*(i+1);}
        return q;
    }

    //###############//
    //### RESULTS ###//
    //###############//

    std::string json_escape(const std::string& s) {
        std::string out;
        for (char c : s) {
            if (c == '"' || c == '\\') {out += '\\'; out += c;}
            else if (c == '\n') {out += "\\n";}
            else if (static_cast<unsigned char>(c) < 0x20) {out += ' ';}
            else {out += c;}
        }
        return out;
    }

    // One line of results.jsonl. The key identifies a measurement, so a restarted run can skip what is already there.
    struct Record {
        std::string suite, program, variant, structure;
        int atoms = 0;
        std::vector<double> samples;
        std::string status = "ok";
        std::map<std::string, std::string> extra;  // already JSON-encoded values

        std::string key() const {return suite + "/" + program + "/" + variant + "/" + structure;}

        std::string json() const {
            std::ostringstream o;
            o << std::setprecision(9);
            o << "{\"key\":\"" << json_escape(key()) << "\",\"suite\":\"" << suite << "\",\"program\":\"" << program
              << "\",\"variant\":\"" << variant << "\",\"structure\":\"" << structure << "\",\"atoms\":" << atoms
              << ",\"status\":\"" << json_escape(status) << "\",\"samples\":[";
            for (size_t i = 0; i < samples.size(); ++i) {o << (i ? "," : "") << samples[i];}
            o << "]";
            for (const auto& [k, v] : extra) {o << ",\"" << k << "\":" << v;}
            o << "}";
            return o.str();
        }
    };

    fs::path results_file() {return config().out / "results.jsonl";}

    // Keys already measured, from a previous (possibly interrupted) run. Failed measurements are not counted, so they
    // are retried.
    std::set<std::string>& done() {
        static std::set<std::string> keys = [] {
            std::set<std::string> k;
            std::ifstream in(results_file());
            for (std::string line; std::getline(in, line);) {
                if (line.find("\"status\":\"ok\"") == std::string::npos && line.find("\"status\":\"over budget\"") == std::string::npos) {continue;}
                auto p = line.find("\"key\":\"");
                if (p == std::string::npos) {continue;}
                p += 7;
                k.insert(line.substr(p, line.find('"', p) - p));
            }
            return k;
        }();
        return keys;
    }

    void write(const Record& r) {
        fs::create_directories(config().out);
        std::ofstream(results_file(), std::ios::app) << r.json() << std::endl;
        if (r.status == "ok" || r.status == "over budget") {done().insert(r.key());}
        std::cout << "  " << std::left << std::setw(60) << r.key() << " ";
        if (r.samples.empty()) {std::cout << r.status << std::endl; return;}
        std::cout << std::fixed << std::setprecision(4) << *std::min_element(r.samples.begin(), r.samples.end()) << " s (min of "
                  << r.samples.size() << ")" << (r.status == "ok" ? "" : " " + r.status) << std::endl;
    }

    bool already_done(const std::string& key) {
        if (!done().contains(key)) {return false;}
        std::cout << "  " << std::left << std::setw(60) << key << " cached" << std::endl;
        return true;
    }

    // One warm-up run, then config().reps timed runs. If the warm-up alone exceeds the budget, only it is kept, and
    // the record is marked so that the caller can stop the series.
    void measure(Record& r, const std::function<void()>& f) {
        using clock = std::chrono::steady_clock;
        auto t0 = clock::now();
        f();
        double warmup = std::chrono::duration<double>(clock::now() - t0).count();
        if (warmup > config().budget) {
            r.samples = {warmup};
            r.status = "over budget";
            return;
        }
        for (int i = 0; i < config().reps; ++i) {
            auto t = clock::now();
            f();
            r.samples.push_back(std::chrono::duration<double>(clock::now() - t).count());
        }
    }

    //#################//
    //### PROCESSES ###//
    //#################//

    // Run argv in cwd with stdin from /dev/null and stdout+stderr to log. Throws if it cannot start or exits non-zero.
    void run_process(const std::vector<std::string>& argv, const fs::path& cwd, const fs::path& log, const std::string& extra_libs = "") {
        std::vector<std::string> env_strings;
        bool ld_set = false;
        for (char** e = environ; *e; ++e) {
            std::string s(*e);
            if (!extra_libs.empty() && s.starts_with("LD_LIBRARY_PATH=")) {
                s = "LD_LIBRARY_PATH=" + extra_libs + ":" + s.substr(16);
                ld_set = true;
            }
            env_strings.push_back(std::move(s));
        }
        if (!extra_libs.empty() && !ld_set) {env_strings.push_back("LD_LIBRARY_PATH=" + extra_libs);}

        std::vector<char*> cargv, cenv;
        for (auto& a : argv) {cargv.push_back(const_cast<char*>(a.c_str()));}
        cargv.push_back(nullptr);
        for (auto& e : env_strings) {cenv.push_back(e.data());}
        cenv.push_back(nullptr);

        // the log is opened after the chdir, so a relative path would be resolved inside cwd
        const std::string abs_cwd = fs::absolute(cwd).string(), abs_log = fs::absolute(log).string();
        posix_spawn_file_actions_t fa;
        posix_spawn_file_actions_init(&fa);
        posix_spawn_file_actions_addchdir_np(&fa, abs_cwd.c_str());
        posix_spawn_file_actions_addopen(&fa, 0, "/dev/null", O_RDONLY, 0);
        posix_spawn_file_actions_addopen(&fa, 1, abs_log.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        posix_spawn_file_actions_adddup2(&fa, 1, 2);

        pid_t pid;
        int rc = posix_spawnp(&pid, cargv[0], &fa, nullptr, cargv.data(), cenv.data());
        posix_spawn_file_actions_destroy(&fa);
        if (rc != 0) {throw std::runtime_error("could not start " + argv[0] + ": " + std::strerror(rc));}
        int status = 0;
        waitpid(pid, &status, 0);
        if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            throw std::runtime_error(argv[0] + " failed (exit " + std::to_string(WIFEXITED(status) ? WEXITSTATUS(status) : -1) + "); see " + log.string());
        }
    }

    // As run_process, but returns the exit code instead of throwing on a non-zero one.
    int run_process_rc(const std::vector<std::string>& argv, const fs::path& cwd, const fs::path& log) {
        try {run_process(argv, cwd, log); return 0;}
        catch (const std::runtime_error& e) {
            std::string msg = e.what();
            auto p = msg.find("(exit ");
            if (p == std::string::npos) {throw;}
            return std::stoi(msg.substr(p + 6));
        }
    }

    std::string read_file(const fs::path& p) {
        std::ifstream in(p);
        return {std::istreambuf_iterator<char>(in), {}};
    }

    // The first line of a program's output that is not decoration, e.g. its version banner.
    std::string first_line_of(const std::vector<std::string>& argv) {
        auto dir = config().out / "work";
        fs::create_directories(dir);
        auto log = dir / "version.log";
        try {run_process(argv, dir, log);} catch (...) {}
        std::istringstream ss(read_file(log));
        for (std::string line; std::getline(ss, line);) {
            if (line.find_first_not_of(" *-\t\r") != std::string::npos) {return line;}
        }
        return "";
    }

    //##############//
    //### INPUTS ###//
    //##############//

    struct SasInput {
        std::string label;
        fs::path pdb, dat;
        bool nm;  // q of the data file is in nm⁻¹
    };

    // manifest.tsv: label, structure path, data path, q unit (A or nm); paths relative to the data root.
    std::vector<SasInput> sas_inputs() {
        std::vector<SasInput> out;
        std::ifstream in(config().manifest);
        if (!in) {return out;}
        for (std::string line; std::getline(in, line);) {
            if (line.empty() || line[0] == '#') {continue;}
            std::stringstream ss(line);
            SasInput s;
            std::string pdb, dat, unit;
            std::getline(ss, s.label, '\t'); std::getline(ss, pdb, '\t'); std::getline(ss, dat, '\t'); std::getline(ss, unit, '\t');
            s.pdb = config().data / pdb;
            s.dat = config().data / dat;
            s.nm = unit == "nm";
            out.push_back(std::move(s));
        }
        return out;
    }

    bool is_water(const std::string& line) {
        static const std::set<std::string> water = {"HOH", "WAT", "SOL", "DOD", "TIP", "TIP3"};
        auto res = line.size() >= 20 ? line.substr(17, 3) : "";
        res.erase(std::remove(res.begin(), res.end(), ' '), res.end());
        return water.contains(res);
    }

    // Hydrogen or deuterium, from the element column or, where that is empty, the atom name without its leading digits.
    bool is_hydrogen(const std::string& line) {
        std::string element = line.size() >= 78 ? line.substr(76, 2) : "";
        element.erase(std::remove(element.begin(), element.end(), ' '), element.end());
        if (!element.empty()) {return element == "H" || element == "D";}
        std::string name = line.size() >= 16 ? line.substr(12, 4) : "";
        auto p = name.find_first_not_of(" 0123456789");
        return p != std::string::npos && (name[p] == 'H' || name[p] == 'D');
    }

    // Write only the heavy, non-solvent ATOM records, so that every program sees the same atoms (they disagree on
    // explicit hydrogens and HETATMs). Returns the number of atoms.
    int strip_structure(const fs::path& src, const fs::path& dst) {
        std::ifstream in(src);
        std::ofstream out(dst);
        int n = 0;
        for (std::string line; std::getline(in, line);) {
            if (!line.empty() && line.back() == '\r') {line.pop_back();}
            if (!line.starts_with("ATOM  ") || is_water(line) || is_hydrogen(line)) {continue;}
            out << line << "\n";
            ++n;
        }
        out << "END\n";
        return n;
    }

    // Copy the data file; one with q in nm⁻¹ is converted to Å⁻¹ and written without its header, so that every program
    // reads the same q values rather than each guessing the unit.
    void prepare_data(const SasInput& s, const fs::path& dst) {
        if (!s.nm) {fs::copy_file(s.dat, dst, fs::copy_options::overwrite_existing); return;}
        std::ifstream in(s.dat);
        std::ofstream out(dst);
        out << std::setprecision(10);
        for (std::string line; std::getline(in, line);) {
            std::istringstream ss(line);
            double q, I, err;
            if (!(ss >> q >> I >> err)) {continue;}
            std::string rest;
            if (ss >> rest) {continue;}  // more than three columns: not a data row
            out << q/10 << " " << I << " " << err << "\n";
        }
    }

    // fcc silver sphere of radius R Å, centred on a lattice site.
    int write_particle(double R, const fs::path& path) {
        constexpr double a = 4.0857; // Å, Ag at room temperature
        constexpr std::array<std::array<double, 3>, 4> basis = {{{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}};
        int n = static_cast<int>(std::ceil(R/a)) + 1;
        std::vector<std::array<double, 3>> atoms;
        for (int i = -n; i <= n; ++i) for (int j = -n; j <= n; ++j) for (int k = -n; k <= n; ++k) {
            for (const auto& b : basis) {
                double x = (i + b[0])*a, y = (j + b[1])*a, z = (k + b[2])*a;
                if (x*x + y*y + z*z <= R*R) {atoms.push_back({x, y, z});}
            }
        }
        fs::create_directories(path.parent_path());
        std::ofstream out(path);
        out << atoms.size() << "\nAg fcc sphere, R = " << R << " A, a = " << a << " A\n" << std::fixed << std::setprecision(6);
        for (const auto& p : atoms) {out << "Ag " << p[0] << " " << p[1] << " " << p[2] << "\n";}
        return static_cast<int>(atoms.size());
    }

    std::string particle_label(double R) {
        std::ostringstream o;
        o << "Ag_R" << R;
        return o.str();
    }

    void write_iq(const fs::path& path, const std::vector<double>& q, const std::vector<double>& I) {
        fs::create_directories(path.parent_path());
        std::ofstream out(path);
        out << std::setprecision(17);
        for (size_t i = 0; i < q.size(); ++i) {out << q[i] << " " << I[i] << "\n";}
    }

    //################//
    //### SETTINGS ###//
    //################//

    // Called before every construction. The thread count only takes effect if it is set before the global pool is
    // first used, which this guarantees since every in-process series starts here.
    void quiet_settings() {
        settings::general::threads = config().threads;
        settings::general::verbose = false;
        settings::general::warnings = false;
        settings::general::offline = true;
        settings::molecule::allow_unknown_atoms = true;
        settings::molecule::allow_unknown_residues = true;
    }

    // Whether the GPU backend actually loads. AUSAXS falls back to the CPU without failing, which must not end up
    // labelled as a GPU measurement.
    bool gpu_usable() {
        static bool usable = [] {
            settings::general::gpu = true;
            return gpu::GPULoader::available();
        }();
        return usable;
    }

    void write_meta() {
        fs::create_directories(config().out);
        std::string cpu;
        {
            std::ifstream in("/proc/cpuinfo");
            for (std::string line; std::getline(in, line);) {
                if (line.starts_with("model name")) {cpu = line.substr(line.find(':') + 2); break;}
            }
        }
        std::ofstream out(config().out / "meta.json");
        out << "{\n"
            << " \"ausaxs_commit\": \"" << AUSAXS_BENCH_COMMIT << "\",\n"
            << " \"cpu\": \"" << json_escape(cpu) << "\",\n"
            << " \"hardware_threads\": " << std::thread::hardware_concurrency() << ",\n"
            << " \"threads\": " << config().threads << ",\n"
            << " \"reps\": " << config().reps << ",\n"
            << " \"gpu\": \"" << json_escape(config().gpu && gpu_usable() ? gpu::GPULoader::device_name() : "") << "\",\n"
            << " \"pepsi\": \"" << json_escape(config().pepsi ? first_line_of({config().pepsi.exe, "--version"}) : "") << "\",\n"
            << " \"crysol\": \"" << json_escape(config().crysol.exe) << "\",\n"
            << " \"foxs\": \"" << json_escape(config().foxs.exe) << "\"\n"
            << "}\n";
    }
}

//##################//
//### SAS SERIES ###//
//##################//

// Wall-clock time of a complete fit to experimental data, from process start to written output, for every program on
// every structure of the manifest. All run in their default modes on the same stripped structure and data file.
TEST_CASE("paper: SAS fit, end to end", "[.][paper][sas]") {
    const auto& c = config();
    auto inputs = sas_inputs();
    if (inputs.empty()) {SKIP("no manifest at " + c.manifest.string());}
    write_meta();

    // program -> (argv template, tool); {pdb} and {dat} are replaced by the prepared copies
    struct Program {std::string name; std::vector<std::string> argv; Tool tool; std::string fallback_marker;};
    const std::vector<std::string> ausaxs_common = {c.ausaxs, "fit", "{pdb}", "{dat}", "--offline", "-o", "out/", "--allow-unknown-atoms", "--allow-unknown-residues"};
    std::vector<Program> programs;
    {
        auto cpu = ausaxs_common; cpu.push_back("--no-gpu");
        programs.push_back({"ausaxs-cpu", cpu, {c.ausaxs, ""}, ""});
        if (c.gpu) {
            auto gpu = ausaxs_common; gpu.push_back("--gpu");
            programs.push_back({"ausaxs-gpu", gpu, {c.ausaxs, ""}, "no usable GPU backend"});
        }
    }
    if (c.pepsi)  {programs.push_back({"pepsi-saxs", {c.pepsi.exe, "{pdb}", "{dat}", "-o", "out.fit"}, c.pepsi, ""});}
    if (c.foxs)   {programs.push_back({"foxs", {c.foxs.exe, "{pdb}", "{dat}"}, c.foxs, ""});}
    if (c.crysol) {programs.push_back({"crysol", {c.crysol.exe, "{dat}", "{pdb}", "--constant", "--implicit-hydrogen", "1"}, c.crysol, ""});}

    std::set<std::string> stopped; // programs whose previous structure exceeded the budget
    for (const auto& s : inputs) {
        if (!fs::exists(s.pdb) || !fs::exists(s.dat)) {WARN("missing inputs for " + s.label); continue;}
        auto base = c.out / "work" / "sas" / s.label;
        fs::create_directories(base);
        int atoms = strip_structure(s.pdb, base / (s.label + ".pdb"));
        std::cout << "\n" << s.label << ": " << atoms << " atoms" << std::endl;

        for (const auto& p : programs) {
            Record r{.suite="sas", .program=p.name, .variant="cli", .structure=s.label, .atoms=atoms};
            if (already_done(r.key())) {continue;}
            if (stopped.contains(p.name)) {r.status = "skipped: budget"; write(r); continue;}

            // a private directory per program, since they all scatter output files next to their inputs
            auto cwd = base / p.name;
            fs::create_directories(cwd);
            auto pdb = cwd / (s.label + ".pdb"), dat = cwd / (s.label + ".dat");
            fs::copy_file(base / (s.label + ".pdb"), pdb, fs::copy_options::overwrite_existing);
            prepare_data(s, dat);
            std::vector<std::string> argv;
            for (auto a : p.argv) {
                if (a == "{pdb}") {a = pdb.filename().string();}
                else if (a == "{dat}") {a = dat.filename().string();}
                argv.push_back(a);
            }
            auto log = cwd / "run.log";
            try {
                measure(r, [&] {run_process(argv, cwd, log, p.tool.libs);});
                if (!p.fallback_marker.empty() && read_file(log).find(p.fallback_marker) != std::string::npos) {
                    r.samples.clear();
                    r.status = "skipped: GPU fallback";
                }
            } catch (const std::exception& e) {
                r.samples.clear();
                r.status = std::string("failed: ") + e.what();
            }
            if (r.status == "over budget") {stopped.insert(p.name);}
            write(r);
        }
    }
}

// The same fit inside one process through the library: load the structure, hydrate, build the histograms and fit.
// The difference to the end-to-end series is process start, output and plotting.
TEST_CASE("paper: SAS fit, in-process", "[.][paper][sas]") {
    const auto& c = config();
    auto inputs = sas_inputs();
    if (inputs.empty()) {SKIP("no manifest at " + c.manifest.string());}

    std::vector<bool> variants = {false};
    if (c.gpu && gpu_usable()) {variants.push_back(true);}
    for (const auto& s : inputs) {
        if (!fs::exists(s.pdb) || !fs::exists(s.dat)) {continue;}
        auto dir = c.out / "work" / "sas" / s.label / "ausaxs-lib";
        fs::create_directories(dir);
        auto pdb = dir / (s.label + ".pdb"), dat = dir / (s.label + ".dat");
        int atoms = strip_structure(s.pdb, pdb);
        prepare_data(s, dat);

        for (bool gpu : variants) {
            Record r{.suite="sas", .program=gpu ? "ausaxs-gpu" : "ausaxs-cpu", .variant="library", .structure=s.label, .atoms=atoms};
            if (already_done(r.key())) {continue;}
            quiet_settings();
            settings::general::gpu = gpu;
            double chi2 = 0;
            try {
                measure(r, [&] {
                    data::Molecule molecule(io::File(pdb.string()));
                    molecule.generate_new_hydration();
                    SimpleDataset data(io::ExistingFile(dat.string()));
                    fitter::SmartFitter fitter(data, molecule.get_histogram());
                    auto result = fitter.fit();
                    chi2 = result->fval/result->dof;
                });
                r.extra["chi2"] = std::to_string(chi2);
            } catch (const std::exception& e) {
                r.samples.clear();
                r.status = std::string("failed: ") + e.what();
            }
            write(r);
        }
    }
    settings::general::gpu = false;
}

//#########################//
//### PURE DEBYE SERIES ###//
//#########################//

namespace {
    // Run one measurement of a Python competitor through competitor.py, which times it in its own interpreter and
    // appends the record itself. Returns false once the series should stop (over budget or failed).
    bool run_competitor(const std::string& suite, const std::string& program, const std::string& variant, const std::string& structure,
                        int atoms, const fs::path& xyz, const std::vector<std::string>& args, bool write_iq_file) {
        const auto& c = config();
        std::string key = suite + "/" + program + "/" + variant + "/" + structure;
        if (already_done(key)) {return true;}
        auto dir = c.out / "work" / suite;
        fs::create_directories(dir);
        std::vector<std::string> argv = {
            c.python, c.adapter, "--suite", suite, "--program", program, "--variant", variant, "--structure", structure,
            "--atoms", std::to_string(atoms), "--xyz", fs::absolute(xyz).string(), "--results", fs::absolute(results_file()).string(),
            "--reps", std::to_string(c.reps), "--threads", std::to_string(c.threads), "--budget", std::to_string(c.budget),
            "--n-q", std::to_string(n_q), "--q-step", std::to_string(q_step)
        };
        if (write_iq_file) {argv.push_back("--iq"); argv.push_back(fs::absolute(c.out / "iq" / (program + "_" + variant + "_" + structure + ".dat")).string());}
        argv.insert(argv.end(), args.begin(), args.end());
        auto log = dir / (program + "_" + variant + "_" + structure + ".log");
        int rc = run_process_rc(argv, fs::current_path(), log);
        done().insert(key); // the adapter wrote its own record, whatever the outcome
        std::cout << "  " << std::left << std::setw(60) << key << " " << (rc == 0 ? "done" : rc == 3 ? "over budget" : "failed, see " + log.string()) << std::endl;
        return rc == 0;
    }

    struct Particle {std::string label; fs::path xyz; int atoms;};

    Particle particle(double R) {
        auto label = particle_label(R);
        auto xyz = config().out / "work" / "particles" / (label + ".xyz");
        int atoms = write_particle(R, xyz);
        return {label, xyz, atoms};
    }

    // Load the particle with the vacuum model, i.e. the plain Debye sum with atomic form factors, at bin width bw.
    // Settings are set before construction, as they must be.
    data::Molecule load_vacuum(const fs::path& xyz, double bw, bool gpu) {
        quiet_settings();
        settings::exv::exv_method = settings::exv::ExvMethod::Vacuum;
        settings::axes::bin_width = bw;
        settings::general::gpu = gpu;
        settings::general::threads = config().threads;
        return data::Molecule(io::File(xyz.string()));
    }

    std::vector<double> ausaxs_forward(data::Molecule& molecule, const std::vector<double>& q) {
        molecule.reset_histogram_manager();
        auto I = molecule.get_histogram()->debye_transform(q);
        std::vector<double> out(I.size());
        for (int i = 0; i < static_cast<int>(I.size()); ++i) {out[i] = I.y(i);}
        return out;
    }
}

// Run time of the forward Debye sum and of its gradient against the number of atoms, for fcc silver spheres on the
// shared q grid. Each program is timed from coordinates in memory to I(q); parsing the file is excluded.
TEST_CASE("paper: pure Debye, scaling", "[.][paper][debye]") {
    const auto& c = config();
    write_meta();
    auto q = debye_q();
    std::vector<double> v(q.size(), 1.0); // dL/dI for the gradient: the gradient of sum_q I(q)

    // AUSAXS: forward at the default bin width and a fine one, on CPU and GPU; gradient on CPU
    struct Variant {std::string name; double bw; bool gpu;};
    std::vector<Variant> forward = {{"cpu-bw0.25", 0.25, false}, {"cpu-bw0.05", 0.05, false}};
    if (c.gpu && gpu_usable()) {forward.push_back({"gpu-bw0.25", 0.25, true}); forward.push_back({"gpu-bw0.05", 0.05, true});}

    std::set<std::string> stopped;
    for (double R : c.radii) {
        auto p = particle(R);
        std::cout << "\n" << p.label << ": " << p.atoms << " atoms" << std::endl;

        for (const auto& var : forward) {
            Record r{.suite="debye", .program="ausaxs", .variant="forward-" + var.name, .structure=p.label, .atoms=p.atoms};
            if (already_done(r.key())) {continue;}
            if (stopped.contains(r.variant)) {r.status = "skipped: budget"; write(r); continue;}
            auto molecule = load_vacuum(p.xyz, var.bw, var.gpu);
            measure(r, [&] {ausaxs_forward(molecule, q);});
            if (r.status == "over budget") {stopped.insert(r.variant);}
            write(r);
        }
        {
            Record r{.suite="debye", .program="ausaxs", .variant="gradient-cpu", .structure=p.label, .atoms=p.atoms};
            if (!already_done(r.key()) && !stopped.contains(r.variant)) {
                auto molecule = load_vacuum(p.xyz, 0.25, false);
                measure(r, [&] {hist::debye_raw_vjp(molecule, q, v);});
                if (r.status == "over budget") {stopped.insert(r.variant);}
                write(r);
            }
        }
        settings::general::gpu = false;

        // competitors: default settings apart from biso = 0, on CPU and on the GPU
        std::vector<std::string> devices = {"cpu"};
        if (c.gpu) {devices.push_back(c.device);}
        for (const auto& dev : devices) {
            if (c.competitors.contains("dc")) {
                for (std::string mode : {"forward", "gradient"}) {
                    std::string variant = mode + "-" + dev;
                    if (stopped.contains("dc/" + variant)) {continue;}
                    if (!run_competitor("debye", "debyecalculator", variant, p.label, p.atoms, p.xyz, {"--device", dev, "--mode", mode}, false)) {stopped.insert("dc/" + variant);}
                }
            }
            if (c.competitors.contains("aes")) {
                std::string variant = "forward-" + dev;
                if (stopped.contains("aes/" + variant)) {continue;}
                if (!run_competitor("debye", "aes-debye", variant, p.label, p.atoms, p.xyz, {"--device", dev, "--mode", "forward"}, false)) {stopped.insert("aes/" + variant);}
            }
        }
    }
}

// Accuracy against the exact Debye sum, and accuracy against cost at fixed size. The reference is DebyeCalculator's
// direct (unbinned) sum in double precision, an implementation independent of AUSAXS. Every profile is written to
// iq/; plot_paper.py takes the deviations.
TEST_CASE("paper: pure Debye, accuracy", "[.][paper][accuracy]") {
    const auto& c = config();
    write_meta();
    auto q = debye_q();
    for (double R : c.accuracy_radii) {
        auto p = particle(R);
        std::cout << "\n" << p.label << ": " << p.atoms << " atoms" << std::endl;

        // the reference, timed once: it is O(N² Q) and only exists to be compared against
        if (c.competitors.contains("dc")) {
            run_competitor("accuracy", "debyecalculator", "exact-float64", p.label, p.atoms, p.xyz,
                {"--device", "cpu", "--mode", "forward", "--pair-sum", "direct", "--dtype", "float64", "--reps", "0"}, true);
        }

        for (double bw : c.bin_widths) {
            std::ostringstream name;
            name << "bw" << bw;
            Record r{.suite="accuracy", .program="ausaxs", .variant=name.str(), .structure=p.label, .atoms=p.atoms};
            if (already_done(r.key())) {continue;}
            auto molecule = load_vacuum(p.xyz, bw, false);
            std::vector<double> I;
            measure(r, [&] {I = ausaxs_forward(molecule, q);});
            write_iq(c.out / "iq" / ("ausaxs_" + name.str() + "_" + p.label + ".dat"), q, I);
            write(r);
        }

        if (c.competitors.contains("dc")) {
            // the grid spacing is 2^floor(log2(f/q_max)); f = 0.12 is the float32 default and 0.03 the float64 one
            for (auto [dtype, f] : std::vector<std::pair<std::string, std::string>>{
                {"float32", "0.48"}, {"float32", "0.12"}, {"float64", "0.12"}, {"float64", "0.03"}
            }) {
                run_competitor("accuracy", "debyecalculator", dtype + "-f" + f, p.label, p.atoms, p.xyz,
                    {"--device", "cpu", "--mode", "forward", "--dtype", dtype, "--spread-step-factor", f}, true);
            }
        }
        if (c.competitors.contains("aes")) {
            for (std::string res : {"1.0", "0.5", "0.25"}) {
                run_competitor("accuracy", "aes-debye", "res" + res, p.label, p.atoms, p.xyz,
                    {"--device", "cpu", "--mode", "forward", "--bins-resolution", res}, true);
            }
        }
    }
}
