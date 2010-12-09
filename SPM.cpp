#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * constructor, makes matrix of dimension M
 * @param M dimension of single particle space and dimension of the Matrix
 * @param N Nr of particles
 */
SPM::SPM(int M,int N) : Matrix(M) {

   this->M = M;
   this->N = N;

}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(SPM &spm_copy) : Matrix(spm_copy) {

   this->M = spm_copy.gM();
   this->N = spm_copy.gN();

}

/**
 * destructor
 */
SPM::~SPM(){

}

/**
 * @return nr of particles
 */
int SPM::gN(){

   return N;

}

/**
 * @return dimension of sp space and of matrix
 */
int SPM::gM(){

   return M;

}

ostream &operator<<(ostream &output,SPM &spm_p){

   for(int i = 0;i < spm_p.M;++i)
      for(int j = 0;j < spm_p.M;++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}
