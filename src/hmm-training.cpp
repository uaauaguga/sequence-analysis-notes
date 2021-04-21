#include <iostream>
#include <map>
#include <getopt.h>
#include "hmm.h"
#include "utils.h"

int main(int argc,char * argv[]){
  int opt;
  std::string emission_path = "params/CpG.emission.random.txt";
  std::string transition_path = "params/CpG.transition.random.txt";
  std::string sequence_path;
  static option long_options[] = {
     {"fasta",required_argument,NULL,'f'},
     {"transition",required_argument,NULL,'t'},
     {"emission",required_argument,NULL,'e'},
     {"help",no_argument,NULL,'h'}
  };
  int option_index = 0;
  while((opt=getopt_long(argc,argv,"f:t:e:h",long_options,&option_index))!=-1){
	  switch(opt){
	    case 'h':
        std::cout<<"Usage:bin/hmm-decoding -s [ACGT]+\n";
        return 0;
        break;
	    case 'f':
		    sequence_path = optarg;
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
  if(sequence_path.length()==0){
    std::cout<<"Input fasta path is required ."<<std::endl;
    return 1;
  }
  std::map<std::string,std::string> sequences;
  //emission.print();
  //transition.print();
  printf("Load fasta\n");
  load_fasta(sequence_path,sequences);
  printf("Training ...\n");
  BM_training(sequences,emission,transition);
  return 0;
}