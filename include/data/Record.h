#pragma once

#include <string>
#include <unordered_map>

class Record {
    public: 
        enum class RecordType {HEADER, ATOM, WATER, TERMINATE, FOOTER, NOTYPE};
        virtual ~Record() = default;
        
        virtual void parse_pdb(const std::string& s) = 0;
        virtual RecordType get_type() const = 0;
        virtual std::string as_pdb() const = 0;

        static RecordType get_type(const std::string& s);

    private:
        // Maps PDB types to a Record. Effectively determines how they are treated by the code.
        static const std::unordered_map<std::string, RecordType> type_map;
};