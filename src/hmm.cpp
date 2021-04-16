#include "hmm.h"
#include <iostream>
float** forward(const std::string& observation,
                const Matrix& emission,const Matrix& transition,
                bool normalize = false){
    int L = observation.length();
    int n_states = emission.nrows; // Number of states with start and end states excluded
    int n_symbol = emission.ncols;  // Number of eimited symbols
    float ** alpha = new float*[n_states];
    for(int s=0;s<n_states;s++){
      alpha[s] = new float[L];
      std::string state = transition.index[s];
      alpha[s][0] = transition.get("s",state)*emission.get(state,observation.substr(0,1));
    }
    // Forward dynamic programming 
    float scaling_factor;
    for(int p=1;p<L;p++){
      scaling_factor = 0;
      for(int s_ = 0;s_<n_states;s_++){
          alpha[s_][p] = 0;
          for(int _s =0;_s<n_states;_s++){
              alpha[s_][p] += alpha[_s][p-1]*transition.values[_s][s_];
          }
          alpha[s_][p] *= emission.values[s_][emission.cn2ci.find(observation.substr(p,1))->second];
          scaling_factor += alpha[s_][p];
      }
      if(normalize){
        for(int s_ = 0;s_<n_states;s_++){
          alpha[s_][p] /= scaling_factor; // Normalize data to prevent numericial under flow
        }
      }
    }
  return alpha;
}


float** backward(const std::string& observation,
                const Matrix& emission,const Matrix& transition,
                bool normalize = false){
    int L = observation.length();
    int n_states = emission.nrows; // Number of states with start and end states excluded
    int n_symbol = emission.ncols;  // Number of eimited symbols
    float **beta = new float*[n_states];
    for(int s=0;s<n_states;s++){
      beta[s] = new float[L];
      std::string state = transition.index[s];
      beta[s][L-1] = transition.get(state,"e");
    }
    float scaling_factor = 0;
    for(int p=L-2;p>=0;p--){
      for(int _s = 0;_s<n_states;_s++){
          beta[_s][p] = 0;
          for(int s_ = 0;s_<n_states;s_++){
            float po_ = emission.values[s_][emission.cn2ci.find(observation.substr(p+1,1))->second];
            float ps_ = transition.values[_s][s_];
            beta[_s][p] += ps_*po_*beta[s_][p+1];
          }
      scaling_factor += beta[_s][p];
      }
      if(normalize){
        for(int _s = 0;_s<n_states;_s++){
          beta[_s][p] /= scaling_factor;
        }
      }
    }
    return beta;
}


float** forward_backward(const std::string& observation,
                const Matrix& emission,const Matrix& transition){
  int L = observation.length();
  int n_states = emission.nrows; // Number of states with start and end states excluded
  int n_symbol = emission.ncols;  // Number of eimited symbols
  float **gamma = new float*[n_states];
  for(int s=0;s<n_states;s++)gamma[s] = new float[L];
  float ** alpha = forward(observation,emission,transition,false);
  float ** beta =  backward(observation,emission,transition,false);
  for(int p=0;p<L;p++){
    float scaling_factor = 0;
    for(int s = 0;s<n_states;s++){
       gamma[s][p] = alpha[s][p]*beta[s][p];
       scaling_factor += gamma[s][p];
    }
    for(int s = 0;s<n_states;s++){
        gamma[s][p] = gamma[s][p]/scaling_factor;
    }
  }
  return gamma;
}

std::string viterbi(const std::string& observation,const Matrix& emission,const Matrix& transition,bool normalize=true){
    int L = observation.length();
    int n_states = emission.nrows; // Number of states with start and end states excluded
    int n_symbol = emission.ncols;  // Number of eimited symbols
    float **psi = new float*[n_states];
    float **delta = new float*[n_states];
    // Initialization
    for(int s=0;s<n_states;s++){
      delta[s] = new float[L];
      psi[s] = new float[L];
      std::string state = transition.index[s];
      delta[s][0] = transition.get("s",state)*emission.get(state,observation.substr(0,1));
      psi[s][0] = 0;
    }
    float max_d,d;
    int max_s;
    for(int p = 1;p<L;p++){
      float scaling_factor = 0;
      for(int s_=0;s_<n_states;s_++){
          max_d = -1;
          for(int _s=0;_s<n_states;_s++){
            d = transition.values[_s][s_]*delta[_s][p-1];
            if(d>max_d){
                max_d = d;
                max_s = _s;
            }
        }
        delta[s_][p] = max_d * emission.values[s_][emission.cn2ci.find(observation.substr(p,1))->second];
        scaling_factor += delta[s_][p];
        psi[s_][p] = max_s;
      }
      if(normalize){
        for(int s_=0;s_<n_states;s_++)
         delta[s_][p] = delta[s_][p]/scaling_factor;
      }
 
    }
    max_d = delta[0][L-1];
    max_s = 0;
    for(int s=1;s<n_states;s++){
      if(delta[s][L-1]>max_d){
        max_s = s;
        max_d = delta[s][L-1];
      }
    }
    std::string decoded = emission.index[max_s];
    int p = L-2;
    while(p>=0){
      max_s = psi[max_s][p+1];
      decoded = emission.index[max_s] + decoded;
      p--;
    }
   return decoded;
  }