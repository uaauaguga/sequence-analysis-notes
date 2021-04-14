#include <iostream>
#include <map>
#include<getopt.h>
#include "utils.h"


float**** alloc_dpmtx(int L_a,int L_b){
  float **** dpmtx; //la*la*lb*lb
  int i,j,k;
  dpmtx = new float***[L_a];
  for(i=0;i<L_a;++i){
    dpmtx[i] = new float**[L_a];      
    for(j=0;j<L_a;++j){
        dpmtx[i][j] = new float*[L_b];
        for(k=0;k<L_b;k++){
            dpmtx[i][j][k] = new float[L_b];
            std::fill_n(dpmtx[i][j][k],L_b,-1); 
            }
      }  
  }
  return dpmtx;
}

void release_dpmtx(float **** dpmtx, int L_a,int L_b){
  int i,j,k;
  for(i=0;i<L_a;++i){
    for(j=0;j<L_a;++j){  
      for(k=0;k<L_b;k++)
        delete dpmtx[i][j][k];
      delete dpmtx[i][j] ;
      }  
   
    delete dpmtx[i];    
  }
  
  delete dpmtx;
}



void init_dpmtx(float **** dpmtx, int L_a,int L_b,int span,float gamma=1){
    int i,j,k,l,d_a,d_b;
    /*
    i: sequence 1, left base
    j: sequence 2, right base
    k: sequence 1, left base
    l: sequence 2, right base
    d_a: sequence 1, distance between left base and right base
    d_b: sequence 2, distance between left base and right base
    */
    for(i=0;i<L_a;++i){
        for(k=0;k<L_b;++k){
            for(d_a=0;d_a<span;d_a++){
                for(d_b=0;d_b<span;d_b++){
                    j = i + d_a;
                    l = k + d_b;
                    dpmtx[i][j][k][l] = abs(d_a-d_b)*gamma;
                }
            }
        }
    }
}


// Refer to https://www.tbi.univie.ac.at/RNA/PMcomp/
void align(std::map<coord,float> bppm_a,std::map<coord,float> bppm_b,
std::string sequence_a,std::string sequence_b,int span,float gamma){
  int L_a = sequence_a.length();
  int L_b = sequence_b.length();
  float**** dpmtx = alloc_dpmtx(L_a,L_b);
  // dpmtx[i][j][k][l]: score of subsequence i..j aligned to subsequence k..l
  float**** dpmtxM = alloc_dpmtx(L_a,L_b);
  // dpmtx[i][j][k][l]: score of subsequence i..j aligned to subsequence k..l, assuming (i,j) and (k,l) pairs
  init_dpmtx(dpmtx,L_a,L_b,span);
  init_dpmtx(dpmtx,L_a,L_b,span,1);
  release_dpmtx(dpmtx,L_a,L_b);
  release_dpmtx(dpmtxM,L_a,L_b);
}


int main(int argc,char * argv[]){
  int opt;
  int span = 3;
  std::string path_a,path_b;
  std::string paired_params = "params/score-paired.txt";
  std::string unpaired_params = "params/score-unpaired.txt";
  float gamma = 1; // Gap penalty
  static option long_options[] = {
	  {"dotplot-A",required_argument,NULL,'a'},
    {"dotplot-B",required_argument,NULL,'b'},
    {"paired-score",required_argument,NULL,'p'},
    {"unpaired-score",required_argument,NULL,'u'},
    {"span",required_argument,NULL,'s'},
    {"help",no_argument,NULL,'h'}
  };
  int option_index = 0;
  while((opt=getopt_long(argc,argv,"a:b:h",long_options,&option_index))!=-1){
	  switch(opt){
	    case 'h':
          std::cout<<"Usage: align-dotplot --dotplot-A first.ps --dotplot-B second.ps\n";
          return 0;
          break;
	    case 'a':
		  path_a = optarg;
	      break;
        case 'b':
          path_b = optarg;
          break;
        case 's':
          span = std::atoi(optarg);
          break;
        case 'p':
          paired_params =  optarg;
          break;
        case 'u':
          unpaired_params =  optarg;
          break;
        default:
		  std::cout<<"Unrecognized argument :"<<opt<<std::endl;
    }
  }

  if(path_a.length()==0||path_b.length()==0){
    std::cout<<"Path of two RNAplfold probabilty dotplot (.ps file) are required ."<<std::endl;
    return 1;
  }else if(!file_exists(path_a)){
    std::cout<<"First dotplot file does not exists."<<std::endl;
    return 2;
  }else if(!file_exists(path_b)){
    std::cout<<"Second dotplot file does not exists."<<std::endl;
    return 3;
  }
  std::cout<<"Input path is : "<<std::endl;
  std::cout<<"First dotplot: "<<path_a<<std::endl;
  std::cout<<"Second dotplot: "<<path_b<<std::endl;
  std::string sequence_a,sequence_b;
  std::map<coord,float> bppm_a,bppm_b;

  std::cout<<"Load data ..."<<std::endl;
  load_ps(path_a,bppm_a,sequence_a);
  load_ps(path_b,bppm_b,sequence_b);
  
  std::cout<<"Input RNA sequence is: "<<std::endl;
  std::cout<<"Sequence a: "<<sequence_a<<std::endl;
  std::cout<<"Sequence b: "<<sequence_b<<std::endl;
  std::cout<<"Load scoring matrix ..."<<std::endl;
  
  std::map<conversion,float> scores;
  load_scoring_mtx(paired_params,scores);
  load_scoring_mtx(unpaired_params,scores);
  
  std::string alphabet = "ACGU";
  for(size_t i = 0;i<alphabet.length();i++){
    scores[std::make_pair(alphabet.substr(i,1),"-")] = -1*gamma;
    scores[std::make_pair("-",alphabet.substr(i,1))] = -1*gamma;
  }

  std::cout<<"Scoring matrix:\n";
  print_scoring_mtx(scores);
  
  //align(bppm_a,bppm_b,sequence_a,sequence_b,span,gamma);
return 0;
}