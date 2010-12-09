#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ifstream;
using std::endl;

#include "include.h"

int DPM::counter = 0;

int **DPM::dp2s;
int ***DPM::s2dp;

/**
 * standard constructor: constructs Matrix object of dimension M*(M - 1)*(M - 2)/6 and
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
DPM::DPM(int M,int N) : Matrix(M*(M - 1)*(M - 2)/6) {

   this->N = N;
   this->M = M;
   this->n = M*(M - 1)*(M - 2)/6;

   if(counter == 0){

      //allocatie van s2dp
      s2dp = new int ** [M];

      for(int i = 0;i < M;++i){

         s2dp[i] = new int * [M];

         for(int j = 0;j < M;++j)
            s2dp[i][j] = new int [M];

      }

      //allocatie van dp2s
      dp2s = new int * [n];

      for(int i = 0;i < n;++i)
         dp2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            for(int c = b + 1;c < M;++c){

               s2dp[a][b][c] = teller;

               dp2s[teller][0] = a;
               dp2s[teller][1] = b;
               dp2s[teller][2] = c;

               ++teller;

            }

   }

   ++counter;

}

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)*(M - 2)/6 and copies the content of dpm_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param dpm_c input DPM to be copied
 */
DPM::DPM(DPM &dpm_c) : Matrix(dpm_c){

   this->N = dpm_c.N;
   this->M = dpm_c.M;
   this->n = M*(M - 1)*(M - 2)/6;

   if(counter == 0){

      //allocatie van s2dp
      s2dp = new int ** [M];

      for(int i = 0;i < M;++i){

         s2dp[i] = new int * [M];

         for(int j = 0;j < M;++j)
            s2dp[i][j] = new int [M];

      }

      //allocatie van dp2s
      dp2s = new int * [n];

      for(int i = 0;i < n;++i)
         dp2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            for(int c = b + 1;c < M;++c){

               s2dp[a][b][c] = teller;

               dp2s[teller][0] = a;
               dp2s[teller][1] = b;
               dp2s[teller][2] = c;

               ++teller;

            }

   }

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists dp2s en s2dp will be deleted.
 */
DPM::~DPM(){

   if(counter == 1){

      for(int i = 0;i < M;++i){

         for(int j = 0;j < M;++j)
            delete [] s2dp[i][j];

         delete [] s2dp[i];

      }

      delete [] s2dp;

      for(int i = 0;i < n;++i)
         delete [] dp2s[i];

      delete [] dp2s;

   }

   --counter;

}

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * DPM(a,b,c,d,e,f) = -DPM(b,a,c,d,e,f) = ...\n\n
 * DPM(a,a,c,d,e,f) = 0\n\n
 * @param a first sp index that forms the dp row index i together with b and c
 * @param b second sp index that forms the dp row index i together with a and c
 * @param c third sp index that forms the dp row index i together with a and b
 * @param d first sp index that forms the dp column index j together with e and z
 * @param e second sp index that forms the dp column index j together with d and z
 * @param z third sp index that forms the dp column index j together with d and e
 * @return the number on place DPM(i,j) with the right phase.
 */
double DPM::operator()(int a,int b,int c,int d,int e,int z) const{

   //eerst kijken of er geen indices gelijk zijn:
   if(a == b || a == c || b == c)
      return 0;

   if(d == e || d == z || e == z)
      return 0;

   //dan kijken wel dp index met welke fase moet genomen worden:
   //eerst voor de i
   int i;

   int phase = 1;

   if(a < b){

      if(b < c)
         i = s2dp[a][b][c];
      else if(c < a)
         i = s2dp[c][a][b];
      else{

         i = s2dp[a][c][b];
         phase *= -1;

      }

   }
   else{

      if(a < c){

         i = s2dp[b][a][c];
         phase *= -1;

      }
      else if(c < b){

         i = s2dp[c][b][a];
         phase *= -1;

      }
      else
         i = s2dp[b][c][a];

   }

   //idem voor j maar met d e z
   int j;

   if(d < e){

      if(e < z)
         j = s2dp[d][e][z];
      else if(z < d)
         j = s2dp[z][d][e];
      else{

         j = s2dp[d][z][e];
         phase *= -1;

      }

   }
   else{

      if(d < z){

         j = s2dp[e][d][z];
         phase *= -1;

      }
      else if(z < e){

         j = s2dp[z][e][d];
         phase *= -1;

      }
      else
         j = s2dp[e][z][d];

   }

   return phase*(*this)(i,j);

}

ostream &operator<<(ostream &output,DPM &dpm_p){

   for(int i = 0;i < dpm_p.n;++i)
      for(int j = 0;j < dpm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << dpm_p.dp2s[i][0] << "\t" << dpm_p.dp2s[i][1] << "\t" << dpm_p.dp2s[i][2]

            << "\t" << dpm_p.dp2s[j][0] << "\t" << dpm_p.dp2s[j][1] << "\t" << dpm_p.dp2s[j][2] << "\t" << dpm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return nr of particles
 */
int DPM::gN(){

   return N;

}

/**
 * @return dimension of sp space
 */
int DPM::gM(){

   return M;

}

/**
 * @return dimension of dp space and of Matrix
 */
int DPM::gn(){

   return n;

}

/**
 * The T1-like (generalized T1) map: maps a TPM object (tpm) on a DPM object (*this)
 * @param A term before the tp part of the map
 * @param B term before the np part of the map
 * @param C term before the sp part of the map
 * @param tpm input TPM
 */
void DPM::T(double A,double B,double C,TPM &tpm){

   SPM spm(C,tpm);

   double ward = 2*B*tpm.trace();

   int a,b,c,d,e,z;

   for(int i = 0;i < n;++i){

      a = dp2s[i][0];
      b = dp2s[i][1];
      c = dp2s[i][2];

      for(int j = i;j < n;++j){

         d = dp2s[j][0];
         e = dp2s[j][1];
         z = dp2s[j][2];

         (*this)(i,j) = 0;

         //no particle stuk
         if(i == j)
            (*this)(i,i) = ward;

         if(a == d){

            //tp stuk
            (*this)(i,j) += A*tpm(b,c,e,z);

            //4 sp stukken
            if(c == z)
               (*this)(i,j) -= spm(e,b);

            if(b == z)
               (*this)(i,j) += spm(c,e);

            if(c == e)
               (*this)(i,j) += spm(b,z);

            if(b == e)
               (*this)(i,j) -= spm(c,z);

         }

         if(b == d){

            //tp stuk
            (*this)(i,j) -= A*tpm(a,c,e,z);

            //2 sp stukken
            if(c == z)
               (*this)(i,j) += spm(a,e);

            if(c == e)
               (*this)(i,j) -= spm(a,z);

         }

         if(b == e){

            //tp stuk
            (*this)(i,j) += A*tpm(a,c,d,z);

            //sp stuk
            if(c == z)
               (*this)(i,j) -= spm(a,d);

         }

         //nu enkel nog tp stukken
         if(c == z)
            (*this)(i,j) += A*tpm(a,b,d,e);

         if(b == z)
            (*this)(i,j) -= A*tpm(a,c,d,e);

         if(c == e)
            (*this)(i,j) -= A*tpm(a,b,d,z);

         if(c == d)
            (*this)(i,j) += A*tpm(a,b,e,z);

      }
   }

   //niet vergeten!
   this->symmetrize();

}
/**
 * The T1-map: maps a TPM object (tpm) on a DPM object (*this). Watch out, like with the TPM::T with option = -1, when M = 2N the Q-like map is singular and the
 * inverse map is undefined, so don't use it.
 * @param option == +1 T1_up map, == -1, inverse T1_down map
 * @param tpm input TPM
 */
void DPM::T(int option, TPM &tpm){

   if(option == 1){

      double a = 1.0;
      double b = 1.0/(N*(N - 1.0));
      double c = 1.0/(N - 1.0);

      this->T(a,b,c,tpm);

   }
   else{//option == -1: inverse T1 down

      //eerst zoeken we de inverse T1 down bar matrix
      //dit is gewoon een inverse Q-like afbeelding:
      TPM bar(M,N);

      double a = 1;
      double b = 3.0/(N*(N - 1.0));
      double c = 0.5/(N - 1.0);

      bar.Q(-1,a,b,c,tpm);

      //dan raisen we dit naar de dp ruimte:
      this->hat(bar);

   }

}

/** 
 * The hat function maps a TPM object tpm to a DPM object (*this) so that bar(this) = tpm,
 * The inverse of the TPM::bar function. It is a T1-like map.
 * @param tpm input TPM
 */
void DPM::hat(TPM &tpm){

   double a = 1.0/(M - 4.0);
   double b = 1.0/((M - 4.0)*(M - 3.0)*(M - 2.0));
   double c = 1.0/((M - 4.0)*(M - 3.0));

   this->T(a,b,c,tpm);

}

/**
 * Deduct from (*this) the T1-map of the unit matrix times a constant (scale)\n\n
 * this -= scale* T1(1) \n\n
 * see notes primal_dual.pdf for more information.
 * @param scale the constant
 */
void DPM::min_tunit(double scale){

   double t = (M*(M - 1.0) - 3.0*N*(M - N))/(N*(N - 1.0));

   scale *= t;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= scale;

}

/**
 * fill the DPM from a file with name filename, where the elements are indicated by their sp-indices
 * @param filename Name of the inputfile
 */
void DPM::in_sp(const char *filename){

   ifstream input(filename);

   double value;

   int a,b,c,d,e,z;

   int i,j;

   while(input >> a >> b >> c >> d >> e >> z >> value){

      i = s2dp[a][b][c];
      j = s2dp[d][e][z];

      std::cout << i << "\t" << j << "\t" << value << endl;

      (*this)(i,j) = value;

   }

   this->symmetrize();

}
