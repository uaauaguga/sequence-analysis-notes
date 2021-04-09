#include <iostream>
#include "trie.h"

node::node(){
  for(int i=0;i<4;i++){
      next[i] = NULL;
  }
  end = false;
}


trie::trie(){
    root = new node();
}

trie::trie(std::map<std::string,std::string> &db){
    root = new node();
    for(std::map<std::string,std::string>::iterator it = db.begin();
    it!=db.end();it++){
        this->insert(it->second);
    }
}


void trie::traverse(){
    node * p = root;
    while(p){
      int i = 0;
      while(i<4){
          if(p->next[i])break;
          i++;
        }
        if(i==4)break;
        p = p->next[i];
        if(p)std::cout<<alphabets[i];
    }
    std::cout<<std::endl;
}

void trie::insert(std::string s){
   node * p = root;
   char c;
   bool end = false;
   for(size_t i=0;i<s.length();i++){
        c = s[i];
        if(p->next[lut[c]]){
           p = p->next[lut[c]];
           p->end = false;
        }else{
           p->next[lut[c]] = new node;
           p = p->next[lut[c]];
           if(i==s.length()-1){
               p->end = true;
           }
        }
   }
}


bool trie::query(std::string &s){
   node * p = root;
   char c;
   bool end = false;
   for(size_t i=0;i<s.length();i++){
        c = s[i];
        if(p->next[lut[c]]){
           p = p->next[lut[c]];
        }else{
           return false;
        }
   }
   return true;   
}
