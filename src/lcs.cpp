#include<stdio.h>
#include<stdlib.h>
#include<getopt.h>
#include<string>
#include<iostream>
#include<map>
#include<vector>


// Get gapped longest common subsequence of two string
// Equivalent to an global alignment where 
// match score is set to 1
// Not allow mismatch, gap score is 0
int ** lcs(std::string& s1,std::string& s2){
  int m = s1.length();
  int n = s2.length();
  int** dpmt = new int*[m];
  int i,j;
  for(i = 0; i < m; ++i)
    dpmt[i] = new int[n];
  for(i = 1; i < m; ++i)dpmt[i][0]=0;
  for(i = 1; i < n; ++i)dpmt[0][i]=0;
  for(i = 1; i < m; ++i){
    for(j = 1; j < n; ++j){
      if(s1[i]==s2[j]){
        dpmt[i][j] = dpmt[i-1][j-1] + 1;
      }else{
        if(dpmt[i-1][j]>=dpmt[i][j-1]){
          dpmt[i][j] = dpmt[i-1][j];
        }else{
          dpmt[i][j] = dpmt[i][j-1];
        }
       }
    }
  }
 return dpmt;
}

void print_lcs(int **dpmt,std::string &s1,std::string &s2){
  int i = s1.length()-1;
  int j = s2.length()-1;
  std::vector<int> indices1,indices2;
  while(i>0&&j>0){
    if(dpmt[i][j]==dpmt[i-1][j]){
        i--;
    }else if(dpmt[i][j]==dpmt[i][j-1]){
        j--;
    }else{
        indices1.push_back(i);
        indices2.push_back(j);
        i--;j--;
    }
  }
  int l = indices1.size();
  for(int k=l-1;k>=0;k--){
   std::cout<<s1[indices1[k]];
  }
   std::cout<<std::endl;
  for(int k=l-1;k>=0;k--){
   std::cout<<s2[indices2[k]];
  }
   std::cout<<std::endl;
}


void print_matrix(int ** mt, int nrows, int ncols){
  for(int i=0;i<nrows;i++){
      for(int j=0;j<ncols;j++){
        std::cout<<mt[i][j]<<"  ";
      }
      std::cout<<"\n";
  }

}

int main(int argc,char * argv[]){
  int opt;
  std::string s1,s2;

  static option long_options[] = {
	{"first",required_argument,NULL,'f'},
	{"second", required_argument,NULL,'s'},
    {"help",no_argument,NULL,'h'}
  };
  int option_index = 0;

  while((opt=getopt_long(argc,argv,"f:s:h",long_options,&option_index))!=-1){
	  switch(opt){
	    case 'h':
            printf("Usage: lcs --first UAAUAAA --second AGCGUAAUUAGUC\n");
            return 0;
            break;
	    case 'f':
		s1 = optarg;
	        break;
        case 's':
        s2 = optarg;
            break;
        default:
        return 1;
  }
}
  std::cout<<s1<<std::endl;
  std::cout<<s2<<std::endl;
  int ** dpmt = lcs(s1,s2);
  print_matrix(dpmt, s1.length(), s2.length());
  print_lcs(dpmt, s1, s2);
  return 0;
}
