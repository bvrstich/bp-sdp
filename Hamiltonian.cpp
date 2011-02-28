#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

int ***Hamiltonian::xys_alpha;
int **Hamiltonian::alpha_xys;

int Hamiltonian::L;
int Hamiltonian::M;

/**
 * function that allocates and constructs the lists.
 * @param L_in the dimension of the square lattice.
 */
void Hamiltonian::init(int L_in){

   L = L_in;
   M = L*L*2;

   //allocate
   xys_alpha = new int ** [L];

   for(int x = 0;x < L;++x){

      xys_alpha[x] = new int * [L];

      for(int y = 0;y < L;++y)
         xys_alpha[x][y] = new int [2];

   }

   alpha_xys = new int * [M];

   for(int i = 0;i < M;++i)
      alpha_xys[i] = new int [3];

   //construct:
   int alpha = 0;

   for(int x = 0;x < L;++x)
      for(int y = 0;y < L;++y)
         for(int s = 0;s < 2;++s){

            alpha_xys[alpha][0] = x;
            alpha_xys[alpha][1] = y;
            alpha_xys[alpha][2] = s;

            xys_alpha[x][y][s] = alpha;

            alpha++;

         }

}

/**
 * deallocates the lists.
 */
void Hamiltonian::clear(){

   //delete xys_alpha
   for(int x = 0;x < L;++x){

      for(int y = 0;y < L;++y)
         delete [] xys_alpha[x][y];

      delete [] xys_alpha[x];

   }

   delete [] xys_alpha;

   //delete alpha_xys
   for(int i = 0;i < M;++i)
      delete [] alpha_xys[i];

   delete [] alpha_xys;

}

/**
 * print the list
 */
void Hamiltonian::print(){

   for(int i = 0;i < M;++i)
      std::cout << i << "\t" << alpha_xys[i][0]<< "\t" << alpha_xys[i][1]<< "\t" << alpha_xys[i][2] << std::endl;

}

/**
 * construct the hopping matrix
 * @param T on exit this Matrix will contain the hopping matrix elements
 */
void Hamiltonian::construct_T(Matrix &T){

   T = 0.0;

   int i,j;

   for(int x = 0;x < L;++x)
      for(int y = 0;y < L;++y)
         for(int s = 0;s < 2;++s){

            i = xys_alpha[x][y][s];
            j = xys_alpha[(x + 1)%L][y][s];

            T(i,j) = -1.0;

            j = xys_alpha[x][(y + 1)%L][s];

            T(i,j) = -1.0;

         }

}
