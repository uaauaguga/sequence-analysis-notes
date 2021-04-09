#include <iostream>
#include <stack>
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

// Perform preorder traverse
void trie::traverse(){
    node * p = root;
    int idx;
    std::stack<node*> pstk;
    std::stack<int> idxstk; // Keep track of the character
    idxstk.push(-1);
    pstk.push(p);
    while(pstk.size()>0){
      p = pstk.top();pstk.pop();
      idx = idxstk.top();idxstk.pop();
      if(idx>=0)std::cout<<alphabets[idx];
      if(p->end)std::cout<<std::endl;
      int i = 3;
      while(i>=0){
        if(p->next[i]){
            pstk.push(p->next[i]);
            idxstk.push(i);
        }
        i--;
      }
    }
    std::cout<<std::endl;
}


// Insert an entry into the trie
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
