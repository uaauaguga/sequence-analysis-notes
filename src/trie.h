#ifndef TRIE_H
#define TRIE_H
#include<string>
#include<map>
#include<vector>

class node{
    public:
    node * next[4];
    bool end;
    node();
};

class trie{
  public:
  node * root = NULL;
  std::map<char, int> lut = {{'A',0}, {'C',1}, {'G',2}, {'T',3}};
  std::string alphabets = "ACGT";

  trie();
  trie(std::map<std::string,std::string> &);
  void insert(std::string);
  void traverse();
  bool query(std::string &query);
};



#endif