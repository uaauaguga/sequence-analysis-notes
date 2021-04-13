#ifndef UTILS_H
#define UTILS_H
#include<sys/stat.h>
#include<string>
#include<map>
#include<tuple>
void load_fasta(std::string path,std::map<std::string,std::string>&);
typedef std::pair<int,int> coord;
void load_ps(std::string,std::map<coord,float>&,std::string&);
inline bool file_exists(const std::string& name) {
  //https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
#endif
