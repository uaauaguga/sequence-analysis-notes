#include "hmm.h"
#include <iostream>

void forward(float ** alpha,const std::string& observation,
                const Matrix& emission,const Matrix& transition,
                bool normalize = false){
    int L = observation.length();
    int n_states = emission.nrows; // Number of states with start and end states excluded
    int n_symbol = emission.ncols;  // Number of eimited symbols
    for(int s=0;s<n_states;s++){
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
}


void backward(float ** beta,const std::string& observation,
                const Matrix& emission,const Matrix& transition,
                bool normalize = false){
    int L = observation.length();
    int n_states = emission.nrows; // Number of states with start and end states excluded
    int n_symbol = emission.ncols;  // Number of eimited symbols
    for(int s=0;s<n_states;s++){
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
}


void forward_backward(float ** alpha,float** beta,float** gamma,const std::string& observation,
                const Matrix& emission,const Matrix& transition){
  int L = observation.length();
  int n_states = emission.nrows; // Number of states with start and end states excluded
  int n_symbol = emission.ncols;  // Number of eimited symbols
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
}

std::string viterbi(const std::string& observation,
                   const Matrix& emission,const Matrix& transition,
                   bool normalize=true){
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

void get_ksi(float ** alpha,
             float ** beta,
             float *** ksi,
             const std::string& observation,
             const Matrix& emission,
             const Matrix& transition){
        
  // Calculate ksi matrix
  
  float scaling_factor;
  int n_states = transition.ncols;
  int L = observation.length();
  
  std::string obs_;
  for(int p=0;p<L+1;p++){
      scaling_factor = 0;
      for(int _s=0;_s<n_states;_s++){
        for(int s_=0;s_<n_states;s_++){
          if(p==0){
            obs_ = observation.substr(p,1);
            if(transition.index[_s]=="s"){
              ksi[_s][s_][p] = transition.values[_s][s_]*emission.cn2ci.find(obs_)->second;
            }else{
              ksi[_s][s_][p] = 0;
            }
          }else if(p<L){
            obs_ = observation.substr(p,1);
            ksi[_s][s_][p] = alpha[_s][p]*transition.values[_s][s_]*emission.cn2ci.find(obs_)->second*beta[s_][p];
          }else{
            if(transition.index[s_]=="e"){
              ksi[_s][s_][p] = transition.values[_s][s_];
            }else{
              ksi[_s][s_][p] = 0;
            }
         }
         scaling_factor += ksi[_s][s_][p];
        }  
      }
      for(int _s=0;_s<n_states;_s++){
        for(int s_=0;s_<n_states;s_++){
          if(scaling_factor>0){
            ksi[_s][s_][p] = ksi[_s][s_][p]/scaling_factor;
          }
        }
      }
   }  
}

void BM_training(const std::map<std::string,std::string>& observations,
                 Matrix& emission,Matrix& transition){
  int n_states = emission.nrows; // Number of states with start and end states excluded
  int n_symbol = emission.ncols;  // Number of eimited symbols

  float ** alpha,**beta,**gamma; // The forward matrix, backward matrix and forward - backward matrix
  float ***ksi; // ksi[i][j][t] is P(s(t)_i->s(t+1)_j|observation,model)
  alpha = new float*[n_states];
  beta = new float*[n_states];
  gamma = new float*[n_states];
  ksi = new float **[n_states];
  std::string observation;
  // Start BM - training
  int epoch = 0;
  float ** transition_ = new float*[n_states];
  float ** emission_ = new float*[n_states];
  for(int s = 0;s<n_states;s++){
    transition_[s] = new float[n_states];
    emission_[s] = new float[n_symbol];
    ksi[s] = new float*[n_states];
  }
  while(epoch<50){   
     // Run Baum-Welch Training across different sequences
    for(int _s = 0;_s<n_states;_s++){
      for(int s_ = 0;s_<n_states;s_++){
        transition_[_s][s_] = 0;
        }
      for(int o = 0;o<n_symbol;o++){
        emission_[_s][o] = 0;
      }
     }
    std::cout<<"$$$$$$"<<epoch<<"$$$$$"<<std::endl;
    for(std::map<std::string,std::string>::const_iterator it  = observations.begin();
        it!=observations.end();it++){
      observation = it->second;
      int L = observation.length();
      for(int _s = 0;_s<n_states;_s++){
        alpha[_s] = new float[L];
        beta[_s] = new float[L];
        gamma[_s] = new float[L];
        for(int s_ = 0;s_<n_states;s_++){
          ksi[_s][s_] = new float[L+1];
        }
      }
      // Calculate forward and backward matrix
      forward(alpha,observation,emission,transition,false); // Calculate forward variable
      backward(beta,observation,emission,transition,false);  // Calculate backward variable
      get_ksi(alpha,beta,ksi,observation,emission,transition); // Calculate ksi matrix 
      
      for(int p = 0;p<L;p++){
        for(int _s = 0;_s<n_states;_s++){
          gamma[_s][p] = 0;
          for(int s_ = 0;s_<n_states;s_++){         
             gamma[_s][p] += ksi[_s][s_][p+1];
           }
          }
      }


      for(int p = 0;p<L+1;p++){
        for(int _s = 0;_s<n_states;_s++){
          for(int s_ = 0;s_<n_states;s_++){
            transition_[_s][s_] += ksi[_s][s_][p];
          }
        }
      }


      for(int _s = 0;_s<n_states;_s++){
        delete alpha[_s];delete beta[_s];delete gamma[_s];
        alpha[_s] = beta[_s] = gamma[_s] = NULL;
        for(int s_ = 0;s_<n_states;s_++){
          delete ksi[_s][s_];
          ksi[_s][s_] = NULL;
        }
      }
     }
     
     for(int _s = 0;_s<n_states;_s++){
     float transition_scaler = 0;
     for(int s_ = 0;s_<n_states;s_++){transition_scaler += transition_[_s][s_];}
     for(int s_ = 0;s_<n_states;s_++){
        if(transition_scaler>0){
         transition.values[_s][s_] = transition_[_s][s_]/transition_scaler;
        }
       }
     float emission_scaler = 0;
     for(int o = 0;o<n_symbol;o++){emission_scaler += emission_[_s][o];}
     for(int o = 0;o<n_symbol;o++){
       if(emission_scaler>0){
         emission.values[_s][o] = emission_[_s][o]/emission_scaler;
       }
       }
     }
    transition.print();
    emission.print();
    epoch += 1;
  }
}