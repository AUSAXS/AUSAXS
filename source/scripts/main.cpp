// includes
#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <iostream>

// ROOT
// #include <TStyle.h>
// #include <TROOT.h>
// #include <TCanvas.h>
// #include <TH1.h>

// my own includes
// #include "data/Atom.h"
// #include "Protein.h"
// #include "Tools.h"
// #include "plot_style.h"

// using namespace ROOT;
using std::string, std::cout, std::endl;

struct bad_format : public std::exception {
    bad_format(const char* msg) : msg(msg) {}
    const char* what() const throw() {return msg;}
    const char* msg;
};

class SMILES {
    public: 
        static void parse(string input) {
            input = '(' + input + ")"; // now it's a valid branch
            parse_branch(input);
        }

        // parse the full branch contained within a pair of parantheses
        static void parse_branch(string substring) {
            cout << substring << endl;
            if (substring[0] != '(' || substring[substring.size()-1] != ')') {
                string err = "Input \"" + substring + "\"is not a valid SMILES format.";
                throw bad_format(err.c_str());
            }
            size_t index = 1;
            int prev_bindings = 0; // number of bindings to previous element in chain
            int next_bindings = 0; // number of bindings to next element in chain
            while (index != substring.size()) {
                int free_bindings = 0;

                switch(substring[index]) {
                    case 'C': {
                        free_bindings += 4 - prev_bindings;
                        prev_bindings = 1;
                        break;
                    } case 'O': {
                        free_bindings += 2 - prev_bindings;
                        prev_bindings = 1;
                        break;
                    } case 'N': {
                        free_bindings += 3 - prev_bindings;
                        prev_bindings = 1;
                        break;
                    } case '(': {
                        size_t branch_length = substring.find_first_of(')')-index+1;
                        parse_branch(substring.substr(index, branch_length));
                        index += branch_length; // skip everything inside the branch
                        continue;
                    } case ')': {
                        [[fallthrough]];
                    } case '=': {
                        index++;
                        prev_bindings = 2;
                        continue;
                    } default: {
                        string err = "Unrecognized character \'" + std::to_string(substring[index]) + "\' in SMILES string.";
                        throw bad_format(err.c_str());
                    }
                }

                // determine number of bindings to next element of the chain
                if (substring[index+1] == '(') {
                    if (substring[index+2] == '=') {free_bindings -= 2;} // double binding to new branch
                    else {free_bindings -= 1;}                           // single binding to new branch
                } else if (substring[index+1] != ')') {                  // check for end of chain
                    free_bindings -= 1;                                  // single binding to next element in chain
                }

                // insert into the map
                string key = substring[index] + std::to_string(result.size()+1); // strict ascending order
                result.insert({key, free_bindings});
                index++;                
            }
        }

        static string reverse(string input) {
            string output = input;
            size_t n = input.size();
            input += ' ';
            for (int i = 0; i < n; i++) {
                output[i] = input[n-i-1];
            }
            input = output;
            for (int i = 0; i < n; i++) {
                if (output[i] == ')') {
                    int branch_len;
                    for (int j = i; j < n; j++) {
                        if (output[j] == '(') {branch_len = j+1;} // +1 to grab the first element after
                    }
                    for (int j = i; j < i+branch_len; j++) {
                        output[j] = input[branch_len-j+1];
                    }
                    i += branch_len;
                }
            }
            return output;
        }

        static void print() {
            for (const auto& pair : result) {
                std::cout << pair.first << " = " << pair.second << std::endl;
            }
        }

    private:
        static inline std::map<string, int> result = {};
};

int main(int argc, char const *argv[]) {
// VOLUME PDB FILE
    // setting::grid::psc = setting::grid::RadialStrategy;
    // setting::grid::ra = 3;
    // setting::grid::rh = 3;

    // Protein protein(argv[1]);
    // protein.generate_new_hydration();
    // protein.generate_volume_file(argv[2]);

// DATABASE GENERATOR
    // SMILES::parse("CCN(CC)C(=O)COCCOCCNC(=O)COCCOCCNC(=O)CCC(C(=O)O)NC(=O)CCCCCCCCCCCCCCCCCCC(=O)O");
    // cout << "CCCCCCCCCCCCCC(=O)O" << endl;
    // cout << SMILES::reverse("CCCCCCCCCCCCCC(=O)O") << endl;
    // SMILES::parse(SMILES::reverse("CCCCCCCCCCCCCC(=O)O"));
    // SMILES::print();

    string smiles = "CCCCCCCCCCCCCC(=O)O";
    cout << "Parsing " << smiles << endl;
    cout << SMILES::reverse(smiles) << endl;
    SMILES::parse(SMILES::reverse(smiles));
    SMILES::print();
    return 0;
}