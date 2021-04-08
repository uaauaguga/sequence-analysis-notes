#include<iostream>
#include <fstream>
#include "utils.h"

void load_fasta(std::string path,std::map<std::string,std::string>& sequences){
std::ifstream fasta(path);
std::string line,seq_id;
while (std::getline(fasta, line))
{
    if(line[0]=='>'){
      seq_id = line.substr(1);
    }else{
     sequences[seq_id] += line;
    }
}
}
