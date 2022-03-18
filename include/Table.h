#include <vector>
#include <unordered_map>
#include <iostream>

using std::vector;

class Table {
    public: 
        Table() : N(0), M(0), data(0) {}

        Table(unsigned int N, unsigned int M) : N(N), M(M), data(N*M) {}

        const double& index(unsigned int i, unsigned int j) const {return data[M*i + j];}

        double& index(unsigned int i, unsigned int j) {return data[M*i + j];}

        unsigned int N, M;
        vector<double> data;
};

template<typename T, typename Q>
class LookupTable : public Table {
    public: 
        LookupTable() {}

        LookupTable(const vector<T>& rows, const vector<Q>& columns) {
            initialize(rows, columns);
        }

        void initialize(const vector<T>& rows, const vector<Q>& columns) {
            N = rows.size(); M = columns.size();
            data.resize(N*M);

            Tmap.reserve(N);
            for (unsigned int i = 0; i < N; i++) {
                Tmap[rows[i]] = i;
            }

            Qmap.reserve(M);
            for (unsigned int i = 0; i < M; i++) {
                Qmap[columns[i]] = i;
            }
        }

        bool is_empty() const {return data.empty();}

        double lookup_index(int i, int j) const {return index(i, j);}

        double lookup(const T row, const Q col) const {return index(Tmap.at(row), Qmap.at(col));}

        void assign_index(int i, int j, double val) {index(i, j) = val;}

        void assign(const T row, const Q col, double val) {index(Tmap.at(row), Qmap.at(col)) = val;}

    private:
        std::unordered_map<T, unsigned int> Tmap;        
        std::unordered_map<Q, unsigned int> Qmap;        
};