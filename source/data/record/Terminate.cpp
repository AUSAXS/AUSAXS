#ifdef (_MSC_VER)
    #pragma warning(disable:4996) // disable sscanf deprecation warning on MSVC
#endif

#include <data/record/Terminate.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>

#include <string>
#include <iomanip>
#include <sstream>

using namespace data::record;
using std::left, std::right, std::setw;

Terminate::Terminate() = default;

Terminate::Terminate(int serial, const std::string& resName, char chainID, int resSeq, const std::string& iCode) {
    this->serial = serial;
    this->resName = resName;
    this->chainID = chainID;
    this->resSeq = resSeq;
    this->iCode = iCode;
}

Terminate::~Terminate() = default;

RecordType Terminate::get_type() const {return RecordType::TERMINATE;}

std::string Terminate::get() const {return as_pdb();}

void Terminate::parse_pdb(const std::string& s) {
    if (s.size() < 28) {return;} // sometimes the terminate record consists only of "TER   "

    // http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#TER

    //                   RN SE S1 RN S2 CI RS iC
    //                   0     1     2        
    //                   0  6  1  8  0  1  2  6  7  
    const char form[] = "%6c%5c%6c%3c%1c%1c%4c%1c";
    std::string recName = "      ", serial = "     ", space1 = "      ", resName = "   ", space2 = " ", resSeq = "    ", iCode = " ";
    char chainID;
    sscanf(s.c_str(), form, recName.data(), serial.data(), space1.data(), resName.data(), 
        space2.data(), &chainID, resSeq.data(), iCode.data());

    // sanity check
    if (Record::get_type(recName) != RecordType::TERMINATE) {
        throw except::parse_error("Terminate::parse_pdb: input string is not \"TER   \" (" + recName + ").");
    }

    // remove any spaces
    serial = utility::remove_all(serial, " ");
    resSeq = utility::remove_all(resSeq, " ");

    // set all of the properties
    try {
        this->serial = std::stoi(serial);
        this->resName = resName;
        this->chainID = chainID;
        this->resSeq = std::stoi(resSeq);
        this->iCode = iCode;
    } catch (const std::exception&) { // catch conversion errors and output a more meaningful error message
        throw except::parse_error("Terminate::parse_pdb: Invalid field values in line \"" + s + "\".");
    }

    // DEBUG OUTPUT
    // cout << s << endl;
    // cout << this->as_pdb() << endl;
}

std::string Terminate::as_pdb() const {
    std::stringstream ss;
    //                   RN SE S1 RN S2 CI RS iC
    //                   0     1     2        
    //                   0  6  1  8  0  1  2  6  7  
    //           format "%6c%5c%7c%2c%1c%1c%4c%1c"
    ss << left << setw(6) << "TER   "                // 1 - 6
        << right << setw(5) << serial                // 7 - 11
        << "      "                                  // 12 - 17
        << right << setw(3) << resName               // 18 - 20
        << " "                                       // 21
        << left << setw(1) << chainID                // 22
        << right << setw(4) << resSeq                // 23 - 26
        << right << setw(1) << iCode                 // 27
        << setw(53) << "                                                     " // 80
        << std::endl;
    return ss.str();
}

void Terminate::set_serial(const int serial) {this->serial = serial;}

bool Terminate::operator==(const Terminate& rhs) const = default;