#include <stdlib.h>
#include <string>
#include <set>
#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include <getopt.h>
#include "utils.h"
#include "motifs.h"

/* 
 * Random projection for intialization of motif location
 * Calculate projection in different sequence, then merge
 * Follows zeros or one occurence assupmtion
 */
std::map<std::string,motif> random_projection(
  std::map<std::string,std::string>& sequences,
  int width=10,int size=5){
  std::set<int> positions; // size position to selected from width 
  while(positions.size() < size){
      positions.insert(rand()%width);
  }
  std::map<std::string,motif> key_ps_dict,key_dict;
  std::string key;
  for(std::map<std::string,std::string>::iterator sit = sequences.begin();
     sit!=sequences.end();sit++){
         key_dict.clear();
         for(int b = 0;b<sit->second.size()-width;b++){
           key = "";
           for(std::set<int>::iterator pit = positions.begin();
             pit != positions.end();pit++){
             key += sit->second[b+*pit];
            }
          location loc(sit->first,b);
          key_dict[key].locations.push_back(loc);
         }
         for(std::map<std::string,motif>::iterator kit = key_dict.begin();
            kit != key_dict.end();kit++){
             int selected = rand()%(kit->second.locations.size());
             key_ps_dict[kit->first].locations.push_back(kit->second.locations[selected]);
             key_ps_dict[kit->first].width = width;
         }
   }
  return key_ps_dict;
}

// Only keep motifs with top n occurence
void filter_motifs(std::map<std::string,motif>& ms,int n=4){
   std::map<std::string,int> counts;
   std::vector<int> occurence;
   std::map<std::string,motif> top_motifs;
   for( std::map<std::string,motif>::iterator it = ms.begin();
     it != ms.end();it++){
         occurence.push_back(it->second.locations.size());
     }
   std::sort(occurence.begin(), occurence.end(), std::greater<int>());
   int min_count = occurence[occurence.size()>n?n:occurence.size()];

   std::vector<std::string> removed;
   for( std::map<std::string,motif>::iterator it = ms.begin();
     it != ms.end();it++){
         if(it->second.locations.size()<min_count){
             removed.push_back(it->first);
         }
     } 
    for(std::vector<std::string>::iterator it = removed.begin();
       it!=removed.end();it++){
           ms.erase(*it);
       }
}

int main(int argc,char * argv[]){
  int opt;
  std::string fasta = "";
  int width=10;
  int size = 5;

  static option long_options[] = {
	{"fasta",required_argument,NULL,'f'},
	{"width", required_argument,NULL,'w'},
    {"size", required_argument,NULL,'s'},
    {"help",no_argument,NULL,'h'}
  };
  int option_index = 0;

  while((opt=getopt_long(argc,argv,"f:w:s:h",long_options,&option_index))!=-1){
	  switch(opt){
	    case 'h':
            std::cout<<"Usage: em-motif --fasta [fasta file] --width [motif width]"<<std::endl;
            return 0;
            break;
	    case 'f':
		    fasta = optarg;
	        break;
        case 'w':
            width = atoi(optarg);
            break;
        case 's':
            size = atoi(optarg);
            break;
        default:
        return 1;
    }
  }
  if(fasta.size()==0){
    std::cout<<"Sequence input --fasta is required"<<std::endl;
    return 2;
  }
  std::cout<<"Input fasta path: "<<fasta<<std::endl;
  std::cout<<"Motif width: "<<width<<std::endl;
  std::cout<<"Seed size: "<<size<<std::endl;
  std::map<std::string,std::string> sequences;
  load_fasta(fasta,sequences);
  std::map<std::string,motif> motifs = random_projection(sequences,width=width,size=size);
  filter_motifs(motifs);
  for(std::map<std::string,motif>::iterator it = motifs.begin();
  it!=motifs.end();it++){
    std::cout<<it->first<<std::endl;
    it->second.print_sequences(sequences);
    it->second.get_PWM(sequences,"ACGT");
  }
  
  return 0;
}
