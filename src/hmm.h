#ifndef HMM_H
#define HMM_H
#include<string>
#include<vector>
#include "utils.h"
float** forward(const std::string& observation,const Matrix& emmision,const Matrix& transition,bool normalize);
float** backward(const std::string& observation,const Matrix& emmision,const Matrix& transition,bool normalize);
float** forward_backward(const std::string& observation,const Matrix& emmision,const Matrix& transition);
std::string viterbi(const std::string& observation,const Matrix& emmision,const Matrix& transition,bool normalize);
#endif