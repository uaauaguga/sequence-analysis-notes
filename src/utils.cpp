#include<iostream>
#include <fstream>
#include<sstream>
#include "utils.h"

void load_fasta(std::string path,std::map<std::string,std::string>& sequences){
  std::ifstream fasta(path);
  std::string line,seq_id;
  while (std::getline(fasta, line))
  {
    if(line[0]=='>'){
      seq_id = line.substr(1);
    }else{
     sequences[seq_id] += line;
    }
  }
}


std::vector<std::string> split(const std::string& s,char sep)
{
	std::string token;
	std::stringstream ss(s);
	std::vector<std::string> tokens;
	while(std::getline(ss,token,sep))
		tokens.push_back(token);
	return tokens;
}

void load_ps(std::string path,std::map<coord,float>& bppm,std::string& sequence){
  std::ifstream ps(path);
  std::string line;

  // Get RNA sequence
  while (std::getline(ps, line)){
    if(line.substr(0,9)=="/sequence"){
      while(std::getline(ps, line)&&line[0]!=')'){
        if(line[line.length()-1]=='\\'){
            line = line.substr(0,line.length()-1);
        }
       sequence += line;
      }
      break;
    }
  } 

  // Get pairing probability
  std::string bppm_start = "%start of base pair probability data";
  std::vector<std::string> fields;
  std::vector<float> prob_unpaired(sequence.length(),1);
  int i,j;
  float prob;
  while (std::getline(ps, line)){
    if(line==bppm_start){
      while(std::getline(ps, line)&&line!="showpage"){
        fields = split(line,' ');
        i = std::atoi(fields[0].c_str())-1;
        j = std::atoi(fields[1].c_str())-1;
        prob = std::atof(fields[2].c_str());
        bppm[std::make_pair(i,j)] = prob;
        prob_unpaired[i] -= prob;
        prob_unpaired[j] -= prob;
      }
    }
  }

  // Get unpaired probability
  for(int i =0;i<prob_unpaired.size();++i){
    if(prob_unpaired[i]>=0){
      bppm[std::make_pair(i,i)] = prob_unpaired[i];
    }else{
      float psum = 0;
      bppm[std::make_pair(i,i)] = 0;
      std::vector<int> indices;
      for(int j=0;j<prob_unpaired.size();j++){
          auto it =  bppm.find(std::make_pair(i,j));
          if(it!=bppm.end()){
            psum += it->second;
            indices.push_back(j);
          }
      }
      for(std::vector<int>::iterator idxiter = indices.begin();idxiter!=indices.end();idxiter++)
        bppm[std::make_pair(i,*idxiter)];
    }
  }
}


 void load_scoring_mtx(std::string path,std::map<conversion,float> &scoring_mtx){
  std::string a,b;
  float score;
  std::cout<<path<<std::endl;
  std::ifstream fscore(path);
  if (fscore.is_open())
  {
   while (fscore >> a >> b>>score){
    scoring_mtx[std::make_pair(a,b)] = score;
    if(a!=b)
      scoring_mtx[std::make_pair(b,a)] = score;
   }
    fscore.close();
  }else{
    std::cout << "Unable to open file"; 
  }
}

void print_scoring_mtx(const std::map<conversion,float>& scores){
  for(std::map<conversion,float>::const_iterator it = scores.begin();
  it!=scores.end();it++){
    std::cout<<it->first.first<<"\t"<<it->first.second<<"\t"<<it->second<<std::endl;
  }
}


void Matrix::read(std::string path){
  std::cout<<path<<std::endl;
  std::ifstream fdata(path);
  std::string line;
  if (fdata.is_open())
  {
    std::vector<std::vector<std::string>> data;
    int n_lines = 0;
    int n_fields;
    while (std::getline(fdata, line)){
      std::vector<std::string> fields = split(line,'\t');
      if(n_lines==0){
        n_fields = fields.size();
      }else{
        if(fields.size()!=n_fields)
          std::cout<<"Number of fields is inconsistent ."<<std::endl;
      }
      data.push_back(fields);
      n_lines += 1;
    }
    nrows = n_lines - 1;
    ncols = n_fields - 1;
    values = new float*[nrows];


    // Get index names
    for(int r = 1;r<n_lines;r++){
      index.push_back(data[r][0]);
      rn2ri[data[r][0]] = r-1;
    }
      

    // Get column names
    for(int c = 1;c<n_fields;c++){
      column.push_back(data[0][c]);
      cn2ci[data[0][c]] = c-1;
    }
     

    // Load data
    for(int r = 0;r<nrows;r++){
      values[r] = new float [ncols];
      for(int c = 0;c<ncols;c++)
        values[r][c] = std::atof(data[r+1][c+1].c_str());
    }
    fdata.close();

  }else{
    std::cout << "Unable to open file"; 
  }  
}

void Matrix::write(std::ostream& fout,std::string sep="\t"){
  fout<<""<<sep;
  for(int c=0;c < ncols;c++)
    fout<<column[c]<<sep;
  fout<<std::endl;
  for(int r = 0;r<nrows;r++){
    fout<<index[r]<<sep;
    for(int c=0;c < ncols;c++)
      fout<<values[r][c]<<sep;
    if(r!=nrows-1)
      fout<<std::endl;
  }
  
}

void Matrix::save(std::string path){
   std::ofstream fout;
   fout.open(path);
   this->write(fout);
   fout.close();
}
 
void Matrix::print(){
   this->write(std::cout);
   std::cout<<"\n";
}


float Matrix::get(std::string rname,std::string cname) const {
  return values[rn2ri.find(rname)->second][cn2ci.find(cname)->second];
}
 