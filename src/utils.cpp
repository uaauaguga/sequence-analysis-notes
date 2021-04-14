#include<iostream>
#include <fstream>
#include<sstream>
#include <vector>
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


std::vector<std::string> split(const std::string& str,char sep)
{
	std::string s=str;
	std::string token;
	std::stringstream ss(s);
	std::vector<std::string> tokens;
	while(std::getline(ss,token,sep))
		tokens.push_back(token);
	return tokens;
}

void load_ps(std::string path,std::map<coord,float>& bppm,std::string& sequence){
  std::ifstream ps(path);
  std::string line;

  // Get RNA sequence
  while (std::getline(ps, line)){
    if(line.substr(0,9)=="/sequence"){
      while(std::getline(ps, line)&&line[0]!=')'){
        if(line[line.length()-1]=='\\'){
            line = line.substr(0,line.length()-1);
        }
       sequence += line;
      }
      break;
    }
  } 

  // Get pairing probability
  std::string bppm_start = "%start of base pair probability data";
  std::vector<std::string> fields;
  std::vector<float> prob_unpaired(sequence.length(),1);
  int i,j;
  float prob;
  while (std::getline(ps, line)){
    if(line==bppm_start){
      while(std::getline(ps, line)&&line!="showpage"){
        fields = split(line,' ');
        i = std::atoi(fields[0].c_str())-1;
        j = std::atoi(fields[1].c_str())-1;
        prob = std::atof(fields[2].c_str());
        bppm[std::make_pair(i,j)] = prob;
        prob_unpaired[i] -= prob;
        prob_unpaired[j] -= prob;
      }
    }
  }

  // Get unpaired probability
  for(int i =0;i<prob_unpaired.size();++i){
    if(prob_unpaired[i]>=0){
      bppm[std::make_pair(i,i)] = prob_unpaired[i];
    }else{
      float psum = 0;
      bppm[std::make_pair(i,i)] = 0;
      std::vector<int> indices;
      for(int j=0;j<prob_unpaired.size();j++){
          auto it =  bppm.find(std::make_pair(i,j));
          if(it!=bppm.end()){
            psum += it->second;
            indices.push_back(j);
          }
      }
      for(std::vector<int>::iterator idxiter = indices.begin();idxiter!=indices.end();idxiter++)
        bppm[std::make_pair(i,*idxiter)];
    }
  }
}


 void load_scoring_mtx(std::string path,std::map<conversion,float> &scoring_mtx){
  std::string a,b;
  float score;
  std::cout<<path<<std::endl;
  std::ifstream fscore(path);
  if (fscore.is_open())
  {
   while (fscore >> a >> b>>score){
    scoring_mtx[std::make_pair(a,b)] = score;
    if(a!=b)
      scoring_mtx[std::make_pair(b,a)] = score;
   }
    fscore.close();
  }else{
    std::cout << "Unable to open file"; 
  }
}

void print_scoring_mtx(const std::map<conversion,float>& scores){
  for(std::map<conversion,float>::const_iterator it = scores.begin();
  it!=scores.end();it++){
    std::cout<<it->first.first<<"\t"<<it->first.second<<"\t"<<it->second<<std::endl;
  }
}