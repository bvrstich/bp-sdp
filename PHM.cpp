#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ifstream;
using std::endl;

#include "include.h"

int PHM::counter = 0;

int **PHM::ph2s;
int **PHM::s2ph;

/**
 * standard constructor: constructs Matrix object of dimension M*M and
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
PHM::PHM(int M,int N) : Matrix(M*M) {
   
   this->N = N;
   this->M = M;
   this->n = M*M;

   if(counter == 0){

      //allocatie van s2ph
      s2ph = new int * [M];
      s2ph[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2ph[i] = s2ph[i - 1] + M;

      //allocatie van ph2s
      ph2s = new int * [n];

      for(int i = 0;i < n;++i)
         ph2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            s2ph[a][b] = teller;

            ph2s[teller][0] = a;
            ph2s[teller][1] = b;

            ++teller;

         }

   }

   ++counter;

}

/**
 * copy constructor: constructs Matrix object of dimension M*M and copies the content of phm_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(PHM &phm_c) : Matrix(phm_c){

   this->N = phm_c.N;
   this->M = phm_c.M;
   this->n = M*M;

   if(counter == 0){

      //allocatie van sp2tp
      s2ph = new int * [M];
      s2ph[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2ph[i] = s2ph[i - 1] + M;

      //allocatie van tp2sp
      ph2s = new int * [n];

      for(int i = 0;i < n;++i)
         ph2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            s2ph[a][b] = teller;

            ph2s[teller][0] = a;
            ph2s[teller][1] = b;

            ++teller;

         }

   }

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists ph2s en s2ph will be deleted.
 */
PHM::~PHM(){

   if(counter == 1){

      delete [] s2ph[0];
      delete [] s2ph;

      for(int i = 0;i < n;++i)
         delete [] ph2s[i];

      delete [] ph2s;

   }

   --counter;

}

/**
 * access the elements of the matrix in sp mode, 
 * @param a first sp index that forms the ph row index i together with b
 * @param b second sp index that forms the ph row index i together with a
 * @param c first sp index that forms the ph column index j together with d
 * @param d second sp index that forms the ph column index j together with c
 * @return the number on place PHM(i,j)
 */
double &PHM::operator()(int a,int b,int c,int d){

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,PHM &phm_p){

   for(int i = 0;i < phm_p.n;++i)
      for(int j = 0;j < phm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << phm_p.ph2s[i][0] << "\t" << phm_p.ph2s[i][1]

            << "\t" << phm_p.ph2s[j][0] << "\t" << phm_p.ph2s[j][1] << "\t" << phm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return number of particles
 */
int PHM::gN(){

   return N;

}

/**
 * @return number of single particle oribals
 */
int PHM::gM(){

   return M;

}

/**
 * @return dimension of the particle hole space, which is also the dimension of the matrix
 */
int PHM::gn(){

   return n;

}

/**
 * De G map, maps a TPM object on a PHM object.
 * @param option = 1 G_up map is used, = -1 G^{-1}_down map is used
 * @param tpm input TPM
 */
void PHM::G(int option,TPM &tpm){

   SPM spm(M,N);

   if(option == 1)
      spm.constr(1.0/(N - 1.0),tpm);
   else
      spm.constr(1.0/(M - N + 1.0),tpm);

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j = i;j < n;++j){

         c = ph2s[j][0];
         d = ph2s[j][1];

         (*this)(i,j) = -tpm(a,d,c,b);

         if(b == d)
            (*this)(i,j) += spm(a,c);

      }
   }
   
   //nog schalen met 4 voor G^{-1}, door de G down die eigenlijk een factor 4 te groot is
   if(option == -1)
      this->dscal(0.25);

   this->symmetrize();

}

/**
 * Calculate the skew trace, defined as:\n\n
 * sum_{ab} PHM(a,a,b,b)
 * @return the skew trace
 */
double PHM::skew_trace(){

   double ward = 0.0;

   for(int a = 0;a < M;++a)
      for(int b = 0;b < M;++b)
         ward += (*this)(a,a,b,b);

   return ward;

}

/**
 * Deduct from this the G-map of the unit matrix times a constant (scale)\n\n
 * this -= scale* G(1) \n\n
 * see notes primal_dual.pdf for more information.
 * @param scale the constant
 */
void PHM::min_gunit(double scale){

   for(int a = 0;a < M;++a)
      for(int b = 0;b < M;++b)
         (*this)(a,a,b,b) -= scale;

   double g = (M - N)/(N - 1.0);

   scale *= g;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= scale;

}

/**
 * Map a PPHM (pphm) object onto a PHM (*this) object by tracing one pair of indices (see primal_dual.pdf for more info)
 * @param pphm Input PPHM
 */
void PHM::bar(PPHM &pphm){

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j = i;j < n;++j){

         c = ph2s[j][0];
         d = ph2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < M;++l)
            (*this)(i,j) += pphm(l,a,b,l,c,d);

      }
   }

   this->symmetrize();

}

/**
 * fill the phm from a file with name filename, where the elements are indicated by their sp-indices
 * @param filename name of the inputfile
 */
void PHM::in_sp(const char *filename){

   ifstream input(filename);

   double value;

   int a,b,c,d;

   while(input >> a >> b >> c >> d >> value)
      (*this)(a,b,c,d) = value;

   this->symmetrize();

}
