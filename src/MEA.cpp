#include <iostream>
#include <map>
#include<getopt.h>
#include "utils.h"

void print_bppm(const std::map<coord,float>& bppm){
  for(std::map<coord,float>::const_iterator it = bppm.begin();
  it!=bppm.end();it++){
    std::cout<<it->first.first<<"\t"<<it->first.second<<"\t"<<it->second<<std::endl;
  }
}


std::string MEA(const std::map<coord,float>& bppm,
                const std::string &sequence,
                const float gamma, const int span){
  std::string structure;
  // Initialize dynamic programming matrix
  int L = sequence.length();
  float ** dpmt = new float * [L];
  for(int i=0;i<L;i++){
    dpmt[i] = new float[L]; 
    std::fill_n(dpmt[i],L,-1); 
  }
  
  // Dynamic programming for MEA
  int i,j,k,d;
  
  for(d = 0;d<L;++d){
    for(i=0;i+d<L;i++){
      j = i + d;
      if(d<span){
      // Number of minimal spanning nucleotide of a base pair is d
        dpmt[i][j] = 0;
      }else{
      // Devide to two case
      // 1. j may unpaired 
      float unpaired_score = dpmt[i][j-1] + bppm.find(std::make_pair(j,j))->second;
      // 2. one of [i,j-d) paired with j
      float paired_prob,paired_score;
      float max_paired_score = 0;
      int max_k;
      for(k=i;k<j-span;k++){
          std::map<std::pair<int,int>,float>::const_iterator it = bppm.find(std::make_pair(k,j));
          if(it==bppm.end())continue;
          paired_prob = it->second;
          paired_score = dpmt[i][k-1] + dpmt[k+1][j-1] + paired_prob*gamma;
          if(paired_score>max_paired_score){
            max_paired_score = paired_score;
            max_k = k;
          }
       }
       if(max_paired_score>=unpaired_score){
         std::cout<<j<<"paired"<<std::endl;
         dpmt[i][j] = max_paired_score;
       }else{
         dpmt[i][j] = unpaired_score;
       }
      }
    }
  }

  
  for(i=0;i<L;i++){
    for(j=0;j<L;j++)
      std::cout<<dpmt[i][j]<<" ";
      std::cout<<"\n";
  }

  for(i=0;i<L;i++){
      delete dpmt[i];
  }
  delete dpmt;

  return structure;
}


int main(int argc,char * argv[]){
  int opt;
  std::string path;
  float gamma = 2;
  int span = 2;
  static option long_options[] = {
	  {"dotplot",required_argument,NULL,'d'},
    {"gamma",optional_argument,NULL,'g'},
    {"span",optional_argument,NULL,'s'},
    {"help",no_argument,NULL,'h'}
  };
  int option_index = 0;
  while((opt=getopt_long(argc,argv,"d:g:s:h",long_options,&option_index))!=-1){
	  switch(opt){
	    case 'h':
        std::cout<<"Usage: MEA --dotplot dotplot.ps\n";
        return 0;
        break;
	    case 'd':
		    path = optarg;
	      break;
      case 'g':
        gamma = std::atof(optarg);
        break;
      case 's':
        span = std::atoi(optarg);
      default:
		    std::cout<<"Unrecognized argument :"<<opt<<std::endl;
    }
  }
  std::cout<<"Input path is : "<<path<<std::endl;
  std::string sequence;
  std::map<coord,float> bppm;
  std::cout<<"Load data ..."<<std::endl;
  load_ps(path,bppm,sequence);
  std::cout<<"Input RNA sequence is: "<<std::endl;
  std::cout<<sequence<<std::endl;
  std::cout<<"Provided base pairing parobability matrix is: "<<std::endl;
  print_bppm(bppm);
  MEA(bppm,sequence,gamma,span);
return 0;
}