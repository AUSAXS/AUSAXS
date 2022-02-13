#pragma once

// includes
#include <string>
#include <vector>
#include "data/Record.h"

using std::string, std::vector;

class Footer : Record {
  public: 
    Footer() : contents("") {}
    ~Footer() override {}

    /**
     * @brief Get the RecordType of this object.
     * @return Record::Footer
     */
    RecordType get_type() const override {return FOOTER;}

    /**
     * @brief Parse a .pdb format Footer string. This is equivalent to the add method.
     * @param s the .pdb format Footer string.
     */
    void parse_pdb(const string s) override {add(s);}

    /**
     * @brief Get the .pdb format representation of this Footer. This is equivalent to the get method.
     * @return the .pdb format Footer string. 
     */
    string as_pdb() const override {return get();}

    /**
     * @brief Add a Footer line to the internal storage of this Footer. 
     * @param s the Footer line. 
     */
    void add(const string s) {contents += s + "\n";}

    /**
     * @brief Get the .pdb format representation of this Footer.
     * @return the .pdb format Footer string. 
     */
    string get() const {return contents;};

  private: 
    string contents;
};