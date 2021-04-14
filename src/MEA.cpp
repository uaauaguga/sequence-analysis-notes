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


void backtracking(float ** dpmt,
                 const std::map<coord,float>& bppm,
                 int i,int j,
                 std::string& structure,
                 float gamma,int span){
  while(i<j){
    int k;
    float epsilon = 0.001;
    float prob_unpaired = bppm.find(std::make_pair(j,j))->second;
    if(abs(dpmt[i][j]-dpmt[i][j-1]-prob_unpaired)<epsilon){
      j--;
    }else{
      for(k=i;k<j-span;k++){
        std::map<coord,float>::const_iterator it = bppm.find(std::make_pair(k,j));
        if(it==bppm.end())continue;
        float prob_paired = it->second;
        if(abs(dpmt[i][j]-dpmt[i][k-1]-dpmt[k+1][j-1]-prob_paired*gamma)<epsilon){
          structure[k] = '(';
          structure[j] = ')';
          backtracking(dpmt,bppm,i,k-1,structure,gamma,span);
          backtracking(dpmt,bppm,k+1,j-1,structure,gamma,span);
          break;
        }
      }
      return ;
    }
  }
}

std::string MEA(const std::map<coord,float>& bppm,
                const std::string &sequence,
                const float gamma, const int span){
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
         dpmt[i][j] = max_paired_score;
       }else{
         dpmt[i][j] = unpaired_score;
       }
      }
    }
  }

  /*
  for(i=0;i<L;i++){
    for(j=0;j<L;j++)
      std::cout<<dpmt[i][j]<<" ";
      std::cout<<"\n";
  }*/
  

  std::string structure(L,'.');
  backtracking(dpmt,bppm,0,L-1,structure,gamma,span);

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
    {"gamma",required_argument,NULL,'g'},
    {"span",required_argument,NULL,'s'},
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
        break;
      default:
		    std::cout<<"Unrecognized argument :"<<opt<<std::endl;
    }
  }

  if(path.length()==0){
    std::cout<<"Path of RNAplfold probabilty dotplot (.ps file) is required ."<<std::endl;
    return 1;
  }else if(!file_exists(path)){
    std::cout<<"Input file does not exists."<<std::endl;
    return 2;
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
  std::string structure = MEA(bppm,sequence,gamma,span);
  std::cout<<"The MEA structure is : "<<std::endl;
  std::cout<<structure<<std::endl;
return 0;
}