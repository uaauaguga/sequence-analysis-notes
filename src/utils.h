#ifndef UTILS_H
#define UTILS_H
#include<string>
#include<map>
#include<tuple>
void load_fasta(std::string path,std::map<std::string,std::string>&);
typedef std::pair<int,int> coord;
void load_ps(std::string,std::map<coord,float>&,std::string&);
#endif
