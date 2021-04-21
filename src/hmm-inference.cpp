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
        std::cout<<"Usage: bin/hmm-inference -s [ACGT]+\n";
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
  //emission.print();
  int s;
  float** alpha = new float*[transition.ncols];
  for(s=0;s<transition.ncols;s++)
    alpha[s] = new float[observation.length()];
  forward(alpha,observation,emission, transition,false);

  float likelihood = 0;
  for(s=0;s<emission.nrows;s++){
      likelihood += alpha[s][observation.length()-1]*transition.values[s][transition.nrows-1];
  }
  std::cout<<"Likelihood by backward algorithm is: "<<likelihood<<std::endl;

  float** beta = new float*[transition.ncols];
  for(s=0;s<transition.ncols;s++)
    beta[s] = new float[observation.length()];
  backward(beta,observation,emission, transition,false);
  likelihood = 0;
  int symbol = emission.cn2ci.find(observation.substr(0,1))->second;
  for(s=0;s<emission.nrows;s++){
      likelihood += beta[s][0]*transition.values[0][s]*emission.values[s][symbol];
  }
  std::cout<<"Likelihood by backward algorithm is: "<<likelihood<<std::endl;
  return 0;
}