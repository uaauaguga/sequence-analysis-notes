#include<stdio.h>
#include<stdlib.h>
#include<getopt.h>
#include<string>
#include<iostream>
#include<vector>
#include "utils.h"

// Pairwise comparison of two string in STL
void findall(std::string &q,std::string &s,std::vector<size_t>& positions){
    size_t position = s.find(q, 0);
    while(position != std::string::npos)
    {
     positions.push_back(position);
     position = s.find(q,position+1);
    }  
}



void findall_brutal(std::string &q, std::string &s,std::vector<size_t>& positions){
  size_t i=0,j;
  if(s.length()<q.length())return;
  size_t n = s.length() - q.length();
  while(i<n){
    j = 0;
    while(j<q.length()&&(q[j]==s[i+j]))j++;
    if(j==q.length())
      positions.push_back(i);
    i++;
  }
}

void get_Z(std::string& s,std::vector<size_t>& Z){
size_t L,R,k;

 //L: left most position of the interval [L,R]
 //R: right most position of the interval [L,R]
 //k: position relative to the left most of the interval
 L = R = 0;
 size_t n = s.length();
 for(size_t i=0;i<n;i++){
   if(i>R){
    L = R = i;
    while(R<n&&s[R-L] == s[R])
      R++;
    Z[i] = R - L;
    R--;
    }
    else{
    k = i - L;
    if(Z[k]<R-i+1){ // Use precalculated Z
      Z[i] = Z[k];
      }else{
      L = i;
      while(R<n&&s[R-L] == s[R])R++;
      Z[i] = R - L;
      R--;
    } 
  }
 }
}



void findall_Z(std::string &q, std::string &s,std::vector<size_t>& positions){
 std::string concate = q + s;
 size_t n = concate.length();
 std::vector<size_t> Z(n,0); // The Z vector
 get_Z(concate,Z);
 for(size_t i=q.length();i<n;i++)
   if(Z[i]==q.length())
     positions.push_back(i-q.length());
}


// Build pi table for KMP string comparison
// Refer to https://www.geeksforgeeks.org/kmp-algorithm-for-pattern-searching/
void get_pi(std::string& s,std::vector<size_t>& pi){
 size_t i = 1; //Current position relative to input string
 size_t l = 0; //Current length of prefix
 pi[0] = 0;
 while(i<s.length()){
  if(s[i]==s[l]){
    //If match, increase i and l
    pi[i++] = ++l;
  }else{
    if(l!=0){
      // Shift l backward, do not increase i
      l = pi[l-1];
    }else{
      pi[i++] = 0;
    }
  }
 }
}

void findall_KMP(std::string &q, std::string &s,std::vector<size_t>& positions){
  if(s.length()<q.length())return;
  std::vector<size_t> pi(q.length());
  get_pi(q,pi);
  size_t n = s.length() - q.length();
  size_t i=0,j=0;
  while(i<n){
    if(q[j]==s[i]){i++;j++;}
    if(j==q.length()){
      positions.push_back(i-j);
      j = pi[j-1];
      }else{
      if(i < n && q[j] != s[i]){
        if(j!=0){
          j = pi[j-1];
        }else{
          i++;
        }
      }
    }
  }
}



void match(std::map<std::string,std::string> db,std::map<std::string,std::string> query){
  std::map<std::string,std::string>::iterator qit;
  std::map<std::string,std::string>::iterator dbit;
  std::vector<size_t> positions;
  for(qit=query.begin();qit!=query.end();qit++){
   for(dbit=db.begin();dbit!=db.end();dbit++){
    positions.clear();
    findall_KMP(qit->second,dbit->second,positions);
    if(positions.size()>0)
      for(std::vector<size_t>::iterator it = positions.begin();it!=positions.end();it++)
	 std::cout<<qit->first<<"\t"<<dbit->first<<"\t"<<*it<<"\t"<<*it+qit->second.size()<<std::endl; 
   }
  }
}



int main(int argc,char * argv[]){
  int opt;
  std::string path_query,path_db;

  static option long_options[] = {
	{"sequences",required_argument,NULL,'s'},
	{"query", required_argument,NULL,'q'},
    {"help",no_argument,NULL,'h'}
  };
  int option_index = 0;

  while((opt=getopt_long(argc,argv,"s:q:h",long_options,&option_index))!=-1){
	  switch(opt){
	    case 'h':
            printf("Usage: test --sequences db.fa --query query.fa\n");
            return 0;
            break;
	    case 's':
		path_db = optarg;
	        break;
        case 'q':
            path_query = optarg;
            break;
        default:
		//printf("Unrecognized argument : %c.\n",opt);
        return 1;
  }
}

  std::map<std::string,std::string> db;
  std::map<std::string,std::string> query;
  load_fasta(path_db,db);
  load_fasta(path_query,query);
  match(db,query);
return 0;
}
