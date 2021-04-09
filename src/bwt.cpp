#include<getopt.h>
#include<string>
#include<vector>
#include<iostream>
#include<algorithm>
#include<numeric>
#include<map>

/*
 * Add an EOF to string
 * Circulate and sort the sequence
 * Take last column
*/
std::string bwt_encoder(std::string s){
    s = s + "$";
    std::string encoded;
    std::vector<std::string> svec;
    for(int i=0;i<s.length();i++)
        svec.push_back(s.substr(i) + s.substr(0,i));
    std::vector<size_t> idx(svec.size());
    for(size_t i=0;i<svec.size();i++)idx[i]=i;
    std::stable_sort(idx.begin(),idx.end(),[&svec](size_t i1, size_t i2) {return svec[i1] < svec[i2];});
    for(std::vector<size_t>::iterator it = idx.begin();
    it!=idx.end();it++){
        encoded += svec[*it][s.length()-1];
    }
    return encoded;
}


/*
* Prepend the encoded string as the first column, sort strings
* Repeat l times
* Take the one end with EOF as decoded string
*/
std::string bwt_decoder(std::string encoded){
  std::vector<std::string> svec;
  int l = encoded.length();
  for(size_t i=0;i<l;i++){
    svec.push_back(encoded.substr(i,1));
  }
  int i = l-1;
  while(i>0){
    std::stable_sort(svec.begin(),svec.end());
    for(size_t k=0;k<l;k++)
      svec[k] = encoded.substr(k,1) + svec[k];
    i--;
  }
  std::string s = "";
  for(size_t i=0;i<l;i++){
    if(svec[i][l-1]=='$'){
      s = svec[i];
      break;
    }
  }
  return s;
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
  std::cout<<"Input sequence is  : "<<sequence<<std::endl;
  std::string encoded = bwt_encoder(sequence);
  std::cout<<"Encoded sequence is: "<<encoded<<std::endl;
  std::cout<<"Decoded sequence is: "<<bwt_decoder(encoded)<<std::endl;
return 0;
}