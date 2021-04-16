#include <iostream>
#include "utils.h"
int main(int argc,char * argv[]){
  Matrix mtx;
  mtx.read("params/CpG.transition.txt");
  std::cout<<mtx.get("bg","bg")<<std::endl;
  mtx.print();
  mtx.save("test.txt");
return 0;
}
