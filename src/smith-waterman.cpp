#include<stdlib.h>
#include<getopt.h>
#include<string>
#include<iostream>
#include<limits>
#include<map>
#include<vector>

void align(std::string& s1,std::string& s2,int match=0,int mismatch=2,int gap_open=5,int gap_extension=3){
  int L1=s1.length(),L2=s2.length();
  int** match_scores = new int*[L1+1];
  int ** insert_1_scores = new int*[L1+1];
  int ** insert_2_scores = new int*[L1+1];
  int i,j;
  int m_m,m_i1,m_i2,i1_m,i1_i1,i2_m,i2_i2;
  int m,i1,i2;
  for(i = 0; i < m+1; ++i){
    match_scores[i] = new int[L2+1];
    insert_1_scores[i] = new int[L2+1];
    insert_2_scores[i] = new int[L2+1];
  }

  int i_max=L1,j_max=L2;
  int max_score = 0;

  for(i = 0; i <= L1; ++i){
    for(j = 0; j <= L2; ++j){
       if(i==0||j==0){
         match_scores[i][j] = insert_1_scores[i][j] = insert_2_scores[i][j] = 0;
       }else{ 
           // case 1. s1[i] matched to s2[j]
           // Maximum score can only be produced in this case, as gaps are penalized
           int sub = (s1[i-1]==s2[j-1])? match:(-1*mismatch);
           m_m = match_scores[i-1][j-1] + sub;
           i1_m = insert_1_scores[i-1][j-1] + sub;
           i2_m = insert_2_scores[i-1][j-1] + sub;
           m = m_m>i1_m?m_m:i1_m;m = m>i2_m?m:i2_m;
           m = m>0?m:0;
           match_scores[i][j] = m;
           if(m > max_score){
             max_score = m;
             i_max = i;
             j_max = j;
           }
           //case 2. s1[i] matched to a gap
           m_i1 = match_scores[i-1][j] - gap_open;
           i1_i1 = insert_1_scores[i-1][j] - gap_extension;
           i1 = m_i1>i1_i1?m_i1:i1_i1;
           i1 = i1>0?i1:0;
           insert_1_scores[i][j] = i1;
           
           //case 3. s2[j] matched to a gap
           m_i2 = match_scores[i][j-1] - gap_open;
           i2_i2 = insert_2_scores[i][j-1] - gap_extension;
           i2 = m_i2>i2_i2?m_i2:i2_i2;
           i2 = i2>0?i2:0;
           insert_2_scores[i][j] = i2;
          }
       }
    }
 
  i = i_max;
  j = j_max;
  std::string align_1 = "";
  std::string align_2 = "";
  while(match_scores[i][j]>0){
    m = match_scores[i][j];
    i1 = insert_1_scores[i][j];
    i2 = insert_2_scores[i][j];
    if(m>i1&&m>i2){
      align_1 = s1[i-1] + align_1;
      align_2 = s2[j-1] + align_2;
      i--;j--;
    }else{
      if(i1>i2){
        align_1 = "-" + align_1;
        align_2 = s2[j-1] + align_2;
        i--;
      }else{
        align_1 = s1[i-1] + align_1;
        align_2 = "-" + align_2;
        j--;
      }
    }
  }
  std::cout<<align_1<<std::endl;
  std::cout<<align_2<<std::endl;
}

// smith–waterman algorithm for local alignment of two sequences
int main(int argc,char * argv[]){
  int opt;
  int gap_open = 5; // penalty for gap opening
  int gap_extension = 3; // penalty for gap extension
  int mismatch = 2; // penalty for mismatch
  int match = 2; // match bonus
  std::string s1,s2;

  static option long_options[] = {
	  {"first",required_argument,NULL,'f'},
	  {"second", required_argument,NULL,'s'},
    {"gap-open", required_argument,NULL,'g'},
    {"gap-extension", required_argument,NULL,'e'},
    {"mismatch", required_argument,NULL,'m'},
    {"match", required_argument,NULL,'b'},
    {"help",no_argument,NULL,'h'}
  };
  int option_index = 0;

  while((opt=getopt_long(argc,argv,"f:s:g:e:m:b:h",long_options,&option_index))!=-1){
	  switch(opt){
	    case 'h':
            printf("Usage: needleman-wunsch --first UAAUAAA --second AGCGUAAUUAGUC\n");
            return 0;
            break;
	    case 'f':
		    s1 = optarg;
	      break;
      case 's':
        s2 = optarg;
        break;
      case 'g':
        gap_open = std::atoi(optarg);
        break;
      case 'e':
        gap_extension = std::atoi(optarg);
      case 'm':
        mismatch = std::atoi(optarg);
        break;
      case 'b':
        match = std::atoi(optarg);
      default:
        return 1;
  }
}
  std::cout<<"Input sequences are:"<<std::endl;
  std::cout<<"First sequence:  "<<s1<<std::endl;
  std::cout<<"Second sequence: "<<s2<<std::endl;
  std::cout<<"Specified parameters are:"<<std::endl;
  std::cout<<"Mismatch penalty:      "<<mismatch<<std::endl;
  std::cout<<"Match bonus:           "<<match<<std::endl;
  std::cout<<"Gap open penalty:      "<<gap_open<<std::endl;
  std::cout<<"Gap extension penalty: "<<gap_extension<<std::endl;
  align(s1,s2,match,mismatch,gap_open,gap_extension);
  return 0;
}
