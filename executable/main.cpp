#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <fstream>
#include <iostream>
#include <memory>

#include <utility/Constants.h>

using std::string, std::cout, std::endl, std::vector;

void parse(std::string id) {
    if (id.size() < 9 || id.substr(0, 9) != "InChI=1S/") {
        cout << "Fatal error, wrong prefix" << endl;
        exit(1);
    }

    std::map<int, int> atoms;
    unsigned int index = 9;
    unsigned int counter = 1;

    // iterate through the entire first section of the id
    while(index < id.size()) {
        // determine the symbol
        string symbol = string(1, id[index]);
        if (symbol == "/") {
            cout << "Found end of first section." << endl;
            break;
        }
        while(std::islower(id[index+1])) {
            symbol += id[++index];
        }
        cout << "symbol is " << symbol << endl;
        unsigned int charge = constants::charge::atomic.at(symbol);

        // determine the number following the symbol
        string number;
        while(index < id.size()) {
            char c = id[++index];
            if (!std::isdigit(c)) {
                break;
            }
            number += c;
        }
        cout << "number is " << number << endl;
        // check we found a number
        if (number.empty()) {
            cout << "Fatal error - invalid string format." << endl;
            exit(2);
        }
        // check we're not dealing with the H (should be ignored here)
        if (symbol == "H") {
            cout << "Found H, skipping" << endl;
            continue;
        } 
        // everything is ok so we add to the map
        for (unsigned int i = 0; i < std::stoi(number); i++) {
            atoms[counter++] = charge;        
        }
    }
}

int main(int argc, char const *argv[]) {
    string test = "InChI=1S/C21H19ClN6O/c22-19-13-24-21-20(26-18(14-28(19)21)15-2-1-7-23-12-15)25-16-3-5-17(6-4-16)27-8-10-29-11-9-27/h1-7,12-14H,8-11H2,(H,25,26)";
    parse(test);
    return 0;
}