#pragma once

// includes
#include <string>
#include <vector>
#include "data/Record.h"

using std::string, std::vector;

class Header : Record {
public: 
    Header(){}
    ~Header() override {}

    /**
     * @brief Get the RecordType of this object.
     * @return Record::HEADER
     */
    RecordType get_type() const override {return HEADER;}

    /**
     * @brief Parse a .pdb format header string. This is equivalent to the add method.
     * @param s the .pdb format header string.
     */
    void parse_pdb(const string s) override {add(contents);}

    /**
     * @brief Get the .pdb format representation of this Header. This is equivalent to the get method.
     * @return the .pdb format header string. 
     */
    string as_pdb() const override {return get();}

    /**
     * @brief Add a header line to the internal storage of this Header. 
     * @param s the header line. 
     */
    void add(const string s) {contents += s + "\n";}

    /**
     * @brief Get the .pdb format representation of this Header.
     * @return the .pdb format header string. 
     */
    string get() const {return contents;};

private: 
    string contents;
};