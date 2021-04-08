#ifndef MOTIFS_H
#define MOTIFS_H
#include<string>
#include<vector>
#include<map>
#include<iostream>
struct location{
  std::string sequence_id;
  int start;
  location(std::string sequence_id,int start):
  sequence_id(sequence_id),start(start){}
};

struct motif{
  int width;
  float pseudo_count = 0.1;
  std::vector<location> locations;
  float ** PWM = NULL;
  int ** frequency = NULL;
  std::map<char,int> alphabet_lut;
  void print_sequences(std::map<std::string,std::string>);
  void get_PWM(std::map<std::string,std::string>,std::string);
};

#endif
