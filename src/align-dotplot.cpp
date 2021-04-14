#include <iostream>
#include <map>
#include<getopt.h>
#include <limits>
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



void sankoff_dp(float **** dpmtx, const std::string& sequence_a,const std::string& sequence_b,
    std::map<coord,float>& bppm_a,std::map<coord,float>& bppm_b,
    int span,const std::map<conversion,float>& scores){
    int i,j,k,l,d_a,d_b;
    /*
    i: sequence 1, left base
    j: sequence 2, right base
    k: sequence 1, left base
    l: sequence 2, right base
    d_a: sequence 1, distance between left base and right base
    d_b: sequence 2, distance between left base and right base
    */


    int L_a = sequence_a.length();
    int L_b = sequence_b.length();

    // Perform four dymensional dynamic programming
    for(d_a=0;d_a<L_a;d_a++){
      std::cout<<"# running: "<<d_a<<std::endl;
      for(d_b=0;d_b<L_b;d_b++){
        for(i=0;i+d_a<L_a;++i){
          for(k=0;k+d_b<L_b;++k){
            j = i + d_a;
            l = k + d_b;
            if(d_a==0||d_b==0){
              dpmtx[i][j][k][l] = scores.find(std::make_pair(sequence_a.substr(i,1),sequence_b.substr(k,1)))->second;
            }else{
              std::map<std::string,float> scoring;
              float max_score = -1*std::numeric_limits<float>::infinity(); // Score to add to i,j;k,l of the dynamic programming matrix
              std::string max_term; // Term that added to the matrix
              // Gaps in four possible position
              scoring["a_l_d"]  = dpmtx[i][j][k+1][l] + scores.find(std::make_pair("-",sequence_b.substr(k,1)))->second;
              scoring["b_l_d"]  = dpmtx[i+1][j][k][l] + scores.find(std::make_pair(sequence_a.substr(i,1),"-"))->second;
              scoring["a_r_d"] = dpmtx[i][j][k][l-1] + scores.find(std::make_pair("-",sequence_b.substr(l,1)))->second;
              scoring["b_r_d"] = dpmtx[i][j-1][k][l] + scores.find(std::make_pair(sequence_a.substr(j,1),"-"))->second;
              // Match with no base pairing involved
              scoring["l_m"] = dpmtx[i+1][j][k+1][l] + scores.find(std::make_pair(sequence_a.substr(i,1),sequence_b.substr(k,1)))->second;
              scoring["r_m"] = dpmtx[i][j-1][k][l-1] + scores.find(std::make_pair(sequence_a.substr(j,1),sequence_b.substr(l,1)))->second;
              // Match with common base pairing
              if(d_a>span&&d_b>span){
                std::string pair_a = sequence_a.substr(i,1) + sequence_a.substr(j,1);
                std::string pair_b = sequence_b.substr(k,1) + sequence_b.substr(l,1);
                float prob_a = bppm_a.find(std::make_pair(i,j))->second; // Probability that i paired with j in sequence a
                float prob_b = bppm_b.find(std::make_pair(k,l))->second; // Probability that k paired with l in sequence b
                // (i,j) and (k,l) form base pair
                scoring["l_r_p"] = dpmtx[i+1][j-1][k+1][l-1] + prob_a + prob_b + scores.find(std::make_pair(pair_a,pair_b))->second;
                // (i,m) and (k,n) form base pair (m<j,n<l) 
                scoring["l_r_up"] = -1*std::numeric_limits<float>::infinity();
                for(int m=span;m<j;m++){
                  for(int n=span;n<l;n++){
                    float score = dpmtx[i][m][k][n] + prob_a + prob_b + dpmtx[m+1][j][n+1][l];
                    if(score>scoring["l_r_up"])scoring["l_r_up"]=score;
                  }
                }
              }
              for(std::map<std::string,float>::iterator it  = scoring.begin();
              it!=scoring.end();++it){
                if(it->second>max_score){
                  max_score = it->second;
                  max_term = it->first;
                }
              }
              dpmtx[i][j][k][l] = max_score;
              //std::cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<l<<"\t"<<max_term<<std::endl;
            }
          }
        }
      }
    }
}


// Refer to https://www.tbi.univie.ac.at/RNA/PMcomp/
void align(std::map<coord,float>& bppm_a,std::map<coord,float>& bppm_b,
std::string sequence_a,std::string sequence_b,int span,std::map<conversion,float> scores){
  int L_a = sequence_a.length();
  int L_b = sequence_b.length();
  float**** dpmtx = alloc_dpmtx(L_a,L_b);
  // dpmtx[i][j][k][l]: score of subsequence i..j aligned to subsequence k..l
  float**** dpmtxM = alloc_dpmtx(L_a,L_b);
  // dpmtx[i][j][k][l]: score of subsequence i..j aligned to subsequence k..l, assuming (i,j) and (k,l) pairs
  sankoff_dp(dpmtx,sequence_a,sequence_b,bppm_a,bppm_b,span,scores);
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
  
  align(bppm_a,bppm_b,sequence_a,sequence_b,span,scores);
return 0;
}