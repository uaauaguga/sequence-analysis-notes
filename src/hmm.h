#ifndef HMM_H
#define HMM_H
#include<string>
#include<vector>
#include "utils.h"
void forward(float ** alpha,
             const std::string& observation,
             const Matrix& emission,
             const Matrix& transition,
             bool normalize);
void backward(float ** alpha,
              const std::string& observation,
              const Matrix& emission,
              const Matrix& transition,
              bool normalize);
void forward_backward(float** alpha, 
                      float** beta,
                      float** gamma,
                      const std::string& observation,
                      const Matrix& emission,
                      const Matrix& transition);
std::string viterbi(const std::string& observation,
                    const Matrix& emission,
                    const Matrix& transition,bool normalize);
void get_ksi(float** alpha,
             float** beta,
             float*** ksi,
             const std::string& observation,
             const Matrix& emission,
             const Matrix& transition);
void BM_training(const std::map<std::string,std::string>& observations,
                 Matrix& emission,
                 Matrix& transition);
void viterbi_training(const std::map<std::string,std::string>& observations,
                     Matrix& emission,
                     Matrix& transition);
#endif