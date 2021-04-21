#include <iostream>
#include <map>
#include <getopt.h>
#include "hmm.h"

int main(int argc,char * argv[]){
  int opt;
  std::string emission_path = "params/CpG.emission.txt";
  std::string transition_path = "params/CpG.transition.txt";
  std::string observation;
  static option long_options[] = {
     {"sequence",required_argument,NULL,'s'},
     {"transition",required_argument,NULL,'t'},
     {"emission",required_argument,NULL,'e'},
     {"help",no_argument,NULL,'h'}
  };
  int option_index = 0;
  while((opt=getopt_long(argc,argv,"s:t:e:h",long_options,&option_index))!=-1){
	  switch(opt){
	    case 'h':
        std::cout<<"Usage:bin/hmm-decoding -s [ACGT]+\n";
        return 0;
        break;
	    case 's':
		  observation = optarg;
	      break;
        case 't':
          transition_path = optarg;
          break;
        case 'e':
          emission_path = optarg;
          break;
        default:
		  std::cout<<"Unrecognized argument :"<<opt<<std::endl;
    }
  }
  Matrix emission, transition;
  emission.read(emission_path);
  transition.read(transition_path);
  int s,p;
  int L = observation.length();

  // Local decoding
  float **alpha,**beta,**gamma;
  alpha = new float*[transition.ncols];
  beta = new float*[transition.ncols];
  gamma = new float*[transition.ncols];
  for(int s = 0;s<transition.ncols;s++){
    alpha[s] = new float[L];
    beta[s] = new float[L];
    gamma[s] = new float[L];
  }
  forward(alpha,observation,emission,transition,false);
  backward(beta,observation,emission,transition,false);
  forward_backward(alpha,beta,gamma,observation,emission,transition);

  /*
  for(s=0;s<emission.nrows;s++){
    for(p=0;p<L;p++){
      printf("%0.2f ",gamma[s][p]);
    }
    printf("\n");
  }*/

  std::cout<<"Result of posterior decoding:"<<std::endl;
  std::cout<<observation<<std::endl;
  for(p=0;p<L;p++){
    float maxp=0;
    int maxs;
    for(s=0;s<emission.nrows;s++){
      if(gamma[s][p]>maxp){
        maxp=gamma[s][p];
        maxs = s;
      }
    }
    std::cout<<transition.index[maxs];
  }
  printf("\n");

  // Viterbi decoding
  std::cout<<"Result of viterbi decoding:"<<std::endl;
  std::string decoded = viterbi(observation, emission, transition,true);
  std::cout<<observation<<std::endl;
  std::cout<<decoded<<std::endl;
  return 0;
}