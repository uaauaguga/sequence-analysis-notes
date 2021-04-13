#include<getopt.h>
#include<string>
#include<iostream>
#include<map>

void backtracking(int ** dpmt,int i,int j,std::string& structure){
  if(i<j){
    int k;
    if(dpmt[i][j]==dpmt[i][j-1]){
      j--;
    }else{
      for(k=i;k<j;k++){
        if(dpmt[i][j]==dpmt[i][k-1]+dpmt[k+1][j-1]+1){
          structure[k] = '(';
          structure[j] = ')';
          backtracking(dpmt,i,k-1,structure);
          backtracking(dpmt,k+1,j-1,structure);
          break;
        }
      }
    }
  }
}


std::string fold(std::string s){

  int l = s.length();
  // Initialize dynamic programming matrix
  int ** dpmt = new int * [l];
  for(int i=0;i<l;i++){
    dpmt[i] = new int[l];
    std::fill_n(dpmt[i],l,0); 
  }
  // Init pairing score matrix
  // delta(i,j): 1 if nuc[i],nuc[j] could form basepair, else 0
  std::map<char,int> lut;
  std::string alphabet = "ACGU";

   int ** delta = new int* [alphabet.length()];
  for(int i=0;i<alphabet.length();i++){
    lut[alphabet[i]] = i;
    delta[i] = new int[alphabet.length()];
    std::fill_n(delta[i],alphabet.length(),0); 
  }

  delta[lut['A']][lut['U']] = 1;
  delta[lut['U']][lut['A']] = 1;
  delta[lut['G']][lut['U']] = 1;
  delta[lut['U']][lut['G']] = 1;
  delta[lut['G']][lut['C']] = 1;
  delta[lut['C']][lut['G']] = 1; 
  
  /*
  * The dynamic programming recursion
  * V: pairing score
  * if i < j, either i and j form base pair
  * V(i,j) = V(i+1,j-1) + delta(nuc[i],nuc[j])
  * or i not form base pair (i=k) / pair with position between i,j (i<k<j) 
  * V(i,j) = max_{i<=k<j}{V(i,k)+V(k+1,j)}
  */
  // Init dynamic prohramming matrix
  for(int i=0;i<s.length();i++){
    dpmt[i][i] = 0;
    if(i<s.length()-1)
      dpmt[i+1][i] =dpmt[i+1][i] = 0;
  }
  // Filling the matrix upward the diagnal
  // Note for d upward off diagnal, i,j follows j - i = d
  // d = 1..l-1, i = 0..l-d

  int i,j,k,d,kmax;
  int paired_score,unpaired_score;
  int max_unpaired_score;
  for(d=1;d<l;d++){
    for(i=0;i<l-d;i++){
      j = i+d;
      max_unpaired_score = 0;
      paired_score = dpmt[i+1][j-1] + delta[lut[s[i]]][lut[s[j]]];
      
      for(k=i;k<j;k++){
        unpaired_score = dpmt[i][k] + dpmt[k+1][j];
        if(unpaired_score>max_unpaired_score){
          kmax=k;
          max_unpaired_score = unpaired_score;
        }
        if(paired_score>max_unpaired_score){
          dpmt[i][j] = paired_score;
        }else{
          dpmt[i][j] = max_unpaired_score;
        }
      }
    }
  }

  /*
  for(i=0;i<l;i++){
    for(j=0;j<l;j++)
      std::cout<<dpmt[i][j]<<" ";
      std::cout<<"\n";
  }*/
  
  std::string structure(s.length(),'.');
  backtracking(dpmt,0,s.length()-1,structure);
  for(i=0;i<l;i++){
      delete dpmt[i];
  }
  delete dpmt;
  return structure;
}



int main(int argc,char * argv[]){
  int opt;
  std::string sequence;
  static option long_options[] = {
	  {"sequence",required_argument,NULL,'s'},
    {"help",no_argument,NULL,'h'}
  };
  int option_index = 0;
  while((opt=getopt_long(argc,argv,"s:h",long_options,&option_index))!=-1){
	switch(opt){
	    case 'h':
        std::cout<<"Usage: nussinov-folding --sequence CGAUGCUA\n";
        return 0;
        break;
	    case 's':
		    sequence = optarg;
	      break;
      default:
		    printf("Unrecognized argument : %c.\n",opt);
    }
  }
  std::cout<<"Input sequence is:   "<<sequence<<std::endl;
  std::string structure = fold(sequence);
  std::cout<<"Output structure is: "<<structure<<std::endl;
return 0;
}
