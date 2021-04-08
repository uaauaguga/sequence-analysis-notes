#include "motifs.h"
void motif::print_sequences(std::map<std::string,std::string> sequences){
    for(std::vector<location>::iterator it = locations.begin();
       it!=locations.end();it++){
         std::cout<<sequences[it->sequence_id].substr(it->start,width)<<std::endl;
       }
  }

void motif::get_PWM(std::map<std::string,std::string> sequences,
    std::string alphabets = "ACGT"){
    int alphabets_size = alphabets.length();

    // Initialize PWM matrix
    PWM = new float *[alphabets_size];
    for(int i = 0;i<alphabets_size;i++)
        alphabet_lut[alphabets[i]] = i;
    for(int i=0;i<alphabets_size;i++){
        PWM[i] = new float[width];
        std::fill_n(PWM[i], width, 0);
    }

    // Initialize frequency matrix
    frequency = new int *[alphabets_size];
    for(int i=0;i<alphabets_size;i++){
        frequency[i] = new int[width];
        std::fill_n(frequency[i], width, 0);
    }

    // Count ocurrence of each nuc at each position
    std::string seq;
    for(std::vector<location>::iterator it = locations.begin();
       it!=locations.end();it++){
         seq = sequences[it->sequence_id].substr(it->start,width);
         for(int i = 0;i<seq.length();i++){
            frequency[alphabet_lut[seq[i]]][i]++;
         }
    }

    // Calculate PWM with pseudocount
    float count;
    for(int j=0;j<width;j++){
        count = 0;
        for(int i=0;i<alphabets_size;i++){
         count += float(frequency[i][j]);
        }
        count += pseudo_count;
        for(int i=0;i<alphabets_size;i++){
           PWM[i][j] = frequency[i][j]/count;
        }
    }
    for(int i=0;i<alphabets_size;i++){
        for(int j=0;j<width;j++){
         std::cout<<PWM[i][j]<<"  ";
        }
        std::cout<<std::endl;
    }    
}