#ifndef UTILS_H
#define UTILS_H
#include<sys/stat.h>
#include<string>
#include<map>
#include <vector>
#include <utility>
#include<tuple>

// Load fasta file
void load_fasta(std::string path,std::map<std::string,std::string>&);

// An integer pair
typedef std::pair<int,int> coord;
void load_ps(std::string,std::map<coord,float>&,std::string&);

// An string pair
typedef std::pair<std::string,std::string> conversion;
void load_scoring_mtx(std::string,std::map<conversion,float>&);
void print_scoring_mtx(const std::map<conversion,float>&);


inline bool file_exists(const std::string& name) {
  //https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

struct Matrix{
    int ncols,nrows;
    std::vector<std::string> index;
    std::vector<std::string> column;
    std::map<std::string,int> cn2ci;
    std::map<std::string,int> rn2ri;
    float ** values;
    void read(std::string);
    void write(std::ostream&,std::string);
    void print();
    void save(std::string);
    float get(std::string,std::string) const ;
};


#endif
