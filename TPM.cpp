#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

int TPM::counter = 0;

int **TPM::t2s;
int **TPM::s2t;

/**
 * standard constructor: constructs Matrix object of dimension M*(M - 1)/2 and
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and tp basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
TPM::TPM(int M,int N) : Matrix(M*(M - 1)/2) {

   this->N = N;
   this->M = M;
   this->n = M*(M - 1)/2;

   if(counter == 0){

      //allocatie van sp2tp
      s2t = new int * [M];
      s2t[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2t[i] = s2t[i - 1] + M;

      //allocatie van tp2sp
      t2s = new int * [n];

      for(int i = 0;i < n;++i)
         t2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j){

            s2t[i][j] = teller;

            t2s[teller][0] = i;
            t2s[teller][1] = j;

            ++teller;

         }

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j)
            s2t[j][i] = s2t[i][j];

   }

   ++counter;

}

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpm_c
 * if counter == 0, the lists containing the relationship between sp and tp basis.
 * @param tpm_c object that will be copied into this.
 */
TPM::TPM(TPM &tpm_c) : Matrix(tpm_c){

   this->N = tpm_c.N;
   this->M = tpm_c.M;
   this->n = M*(M - 1)/2;

   if(counter == 0){

      //allocatie van sp2tp
      s2t = new int * [M];
      s2t[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2t[i] = s2t[i - 1] + M;

      //allocatie van tp2sp
      t2s = new int * [n];

      for(int i = 0;i < n;++i)
         t2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j){

            s2t[i][j] = teller;

            t2s[teller][0] = i;
            t2s[teller][1] = j;

            ++teller;

         }

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j)
            s2t[j][i] = s2t[i][j];

   }

   ++counter;

}

/**
 * construct from file: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of the file
 * if counter == 0, the lists containing the relationship between sp and tp basis.
 * @param filename name of the input file
 */
TPM::TPM(const char *filename) : Matrix(filename){

   ifstream input(filename);

   input >> this->n;

   int I,J;
   double value;

   //inefficient way of going to the last line of the file:
   for(int i = 0;i < n;++i)
      for(int j = 0;j < n;++j)
         input >> I >> J >> value;

   input >> this->M >> this->N;

  if(counter == 0){

      //allocatie van sp2tp
      s2t = new int * [M];
      s2t[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2t[i] = s2t[i - 1] + M;

      //allocatie van tp2sp
      t2s = new int * [n];

      for(int i = 0;i < n;++i)
         t2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j){

            s2t[i][j] = teller;

            t2s[teller][0] = i;
            t2s[teller][1] = j;

            ++teller;

         }

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j)
            s2t[j][i] = s2t[i][j];

   }

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists t2s en s2t will be deleted.
 * 
 */
TPM::~TPM(){

   if(counter == 1){

      delete [] s2t[0];
      delete [] s2t;

      for(int i = 0;i < n;++i)
         delete [] t2s[i];

      delete [] t2s;

   }

   --counter;

}

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * TPM(a,b,c,d) = -TPM(b,a,c,d) = -TPM(a,b,d,c) = TPM(b,a,c,d)
 * @param a first sp index that forms the tp row index i together with b
 * @param b second sp index that forms the tp row index i together with a
 * @param c first sp index that forms the tp column index j together with d
 * @param d second sp index that forms the tp column index j together with c
 * @return the number on place TPM(i,j) with the right phase.
 */
double TPM::operator()(int a,int b,int c,int d) const{

   if( (a == b) || (c == d) )
      return 0;
   else{

      int i = s2t[a][b];
      int j = s2t[c][d];

      int phase = 1;

      if(a > b)
         phase *= -1;
      if(c > d)
         phase *= -1;

      return phase*(*this)(i,j);

   }

}

ostream &operator<<(ostream &output,TPM &tpm_p){

   for(int i = 0;i < tpm_p.n;++i)
      for(int j = 0;j < tpm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << tpm_p.t2s[i][0] << "\t" << tpm_p.t2s[i][1]

            << "\t" << tpm_p.t2s[j][0] << "\t" << tpm_p.t2s[j][1] << "\t" << tpm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return number of particles
 */
int TPM::gN(){

   return N;

}

/**
 * @return number of sp orbitals
 */
int TPM::gM(){

   return M;

}

/**
 * @return de dimensie of the tp matrix space
 */
int TPM::gn(){

   return n;

}

/**
 * construct the hubbard hamiltonian with on site repulsion U
 * @param U onsite repulsion term
 * @param option == 0 use periodic boundary conditions, == 1 use no pbc
 */
void TPM::hubbard(int option,double U){

   int a,b,c,d;//sp orbitals

   double ward = 1.0/(N - 1.0);

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0;

         if(option == 0){//pbc

            //eerst hopping
            if( (a == c) && ( ( (b + 2)%M == d ) || ( b == (d + 2)%M ) ) )
               (*this)(i,j) -= ward;

            if( (b == c) && ( ( (a + 2)%M == d ) || ( a == (d + 2)%M ) ) )
               (*this)(i,j) += ward;

            if( (b == d) && ( ( (a + 2)%M == c ) || ( a == (c + 2)%M ) ) )
               (*this)(i,j) -= ward;

         }
         else{//no pbc

            //eerst hopping
            if( (a == c) && ( ( (b + 2) == d ) || ( b == (d + 2) ) ) )
               (*this)(i,j) -= ward;

            if( (b == c) && ( ( (a + 2) == d ) || ( a == (d + 2) ) ) )
               (*this)(i,j) += ward;

            if( (b == d) && ( ( (a + 2) == c ) || ( a == (c + 2) ) ) )
               (*this)(i,j) -= ward;

         }

         //on site interaction
         if( (a % 2) == 0 && (c % 2) == 0 )
            if(a == (b - 1) && c == (d - 1) && a == c)
               (*this)(i,j) += U;

      }

   }

   this->symmetrize();

}

/**
 * The Q map
 * @param option = 1, regular Q map , = -1 inverse Q map
 * @param tpm_d the TPM of which the Q map is taken and saved in this.
 */
void TPM::Q(int option,TPM &tpm_d){

   double a = 1;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

/**
 * Thee Q-like map: see primal-dual.pdf for more info (form: Q(A,B,C)(TPM) )
 * @param option = 1, regular Q-like map , = -1 inverse Q-like map
 * @param A factor in front of the two particle piece of the map
 * @param B factor in front of the no particle piece of the map
 * @param C factor in front of the single particle piece of the map
 * @param tpm_d the TPM of which the Q-like map is taken and saved in this.
 */
void TPM::Q(int option,double A,double B,double C,TPM &tpm_d){

   if(option == -1){

      B = (B*A + B*C*M - 2.0*C*C)/( A * (C*(M - 2.0) -  A) * ( A + B*M*(M - 1.0) - 2.0*C*(M - 1.0) ) );
      C = C/(A*(C*(M - 2.0) - A));
      A = 1.0/A;

   }

   SPM spm(M,N);

   //de trace*2 omdat mijn definitie van trace in berekeningen over alle (alpha,beta) loopt
   double ward = B*tpm_d.trace()*2.0;

   //construct de spm met schaling C
   spm.constr(C,tpm_d);

   for(int i = 0;i < n;++i){

      int a = t2s[i][0];
      int b = t2s[i][1];

      for(int j = i;j < n;++j){

         int c = t2s[j][0];
         int d = t2s[j][1];

         (*this)(i,j) = A*tpm_d(i,j);

         if(i == j)
            (*this)(i,i) += ward;

         if(a == c)
            (*this)(i,j) -= spm(b,d);

         if(b == c)
            (*this)(i,j) += spm(a,d);

         if(b == d)
            (*this)(i,j) -= spm(a,c);

      }
   }

   this->symmetrize();

}
/**
 * initialize this onto the unitmatrix with trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = N*(N - 1.0)/(2.0*n);

   for(int i = 0;i < n;++i){

      (*this)(i,i) = ward;

      for(int j = i + 1;j < n;++j)
         (*this)(i,j) = (*this)(j,i) = 0.0;

   }

}

/**
 * orthogonal projection onto the space of traceless matrices
 */
void TPM::proj_Tr(){

   double ward = (this->trace())/(double)n;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= ward;

}

/**
 * Primal hessian map:\n\n
 * Hb = D_1 b D_1 + D_2 Q(b) D_2 + D_3 G(b) D_3 + D_4 T1(b) D_4 + D_5 T2(b) D5 \n\n
 * with D_1, D_2, D_3 and D_3 the P, Q, G, T1 and T2 blocks of the SUP D. 
 * @param b TPM domain matrix, hessian will act on it and the image will be put in this
 * @param D SUP matrix that defines the structure of the hessian map. (see primal-dual.pdf for more info)
 */
void TPM::H(TPM &b,SUP &D){

   this->L_map(D.tpm(0),b);

#ifdef __Q_CON

   //maak Q(b)
   TPM Qb(M,N);
   Qb.Q(1,b);

   TPM hulp(M,N);

   hulp.L_map(D.tpm(1),Qb);

   Qb.Q(1,hulp);

   *this += Qb;

#endif

#ifdef __G_CON

   //maak G(b)
   PHM Gb(M,N);
   Gb.G(1,b);

   PHM hulpje(M,N);

   hulpje.L_map(D.phm(),Gb);

   hulp.G(1,hulpje);

   *this += hulp;

#endif

#ifdef __T1_CON

   DPM T1b(M,N);
   T1b.T(1,b);

   DPM hulp_T1(M,N);

   hulp_T1.L_map(D.dpm(),T1b);

   hulp.T(1,hulp_T1);

   *this += hulp;

#endif

#ifdef __T2_CON

   PPHM T2b(M,N);
   T2b.T(0,b);

   PPHM hulp_T2(M,N);

   hulp_T2.L_map(D.pphm(),T2b);

   hulp.T(hulp_T2);

   *this += hulp;

#endif

   this->proj_Tr();

}

/**
 * Implementation of a linear conjugate gradient algoritm for the solution of the primal Newton equations\n\n
 * H(*this) =  b\n\n 
 * in which H represents the hessian map.
 * @param b righthandside of the equation
 * @param D SUP matrix that defines the structure of the hessian
 * @return return number of iterations needed to converge to the desired accuracy
 */
int TPM::solve(TPM &b,SUP &D){

   *this = 0;

   //de r initialiseren op b
   TPM r(b);

   double rr = r.ddot(r);
   double rr_old,ward;

   //nog de Hb aanmaken ook, maar niet initialiseren:
   TPM Hb(M,N);

   int cg_iter = 0;

   while(rr > 1.0e-10){

      ++cg_iter;

      Hb.H(b,D);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      b.dscal(rr/rr_old);

      //dan r er bijtellen
      b += r;

   }

   return cg_iter;

}


/**
 * The G-map that maps a PHM object onto a TPM object.
 * @param option = 1, G_down - map is used, = -1 G^{-1}_up - map is used.
 * @param phm input PHM 
 */
void TPM::G(int option,PHM &phm){

   SPM spm(M,N);

   if(option == 1)
      spm.constr(1.0/(N - 1.0),phm);
   else
      spm.constr(1.0/(M - N + 1.0),phm);

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = phm(b,d,c,a) - phm(a,d,c,b) - phm(b,c,d,a) + phm(a,c,d,b);

         if(b == d)
            (*this)(i,j) += spm(a,c);

         if(b == c)
            (*this)(i,j) -= spm(a,d);

         if(a == c)
            (*this)(i,j) += spm(b,d);

      }

   }

   //nog schalen met 4 voor G^{-1}
   if(option == -1)
      this->dscal(0.25);

   this->symmetrize();

}

/**
 * ( Overlapmatrix of the U-basis ) - map, maps a TPM onto a different TPM, this map is actually a Q-like map
 * for which the paramaters a,b and c are calculated in primal_dual.pdf. Since it is a Q-like map the inverse
 * can be taken as well.
 * @param option = 1 direct overlapmatrix-map is used , = -1 inverse overlapmatrix map is used
 * @param tpm_d the input TPM
 */
void TPM::S(int option,TPM &tpm_d){

   double a = 1.0;
   double b = 0.0;
   double c = 0.0;

#ifdef __Q_CON

   a += 1.0;
   b += (4.0*N*N + 2.0*N - 4.0*N*M + M*M - M)/(N*N*(N - 1.0)*(N - 1.0));
   c += (2.0*N - M)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __G_CON

   a += 4.0;
   c += (2.0*N - M - 2.0)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __T1_CON

   a += M - 4.0;
   b += (M*M*M - 6.0*M*M*N -3.0*M*M + 12.0*M*N*N + 12.0*M*N + 2.0*M - 18.0*N*N - 6.0*N*N*N)/( 3.0*N*N*(N - 1.0)*(N - 1.0) );
   c -= (M*M + 2.0*N*N - 4.0*M*N - M + 8.0*N - 4.0)/( 2.0*(N - 1.0)*(N - 1.0) );

#endif

#ifdef __T2_CON
   
   a += 5.0*M - 8.0;
   b += 2.0/(N - 1.0);
   c += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M)/(2.0*(N - 1.0)*(N - 1.0));

#endif

   this->Q(option,a,b,c,tpm_d);

}

/**
 * Deduct the unitmatrix times a constant (scale) from this.\n\n
 * this -= scale* 1
 * @param scale the constant
 */
void TPM::min_unit(double scale){

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= scale;

}

/**
 * Deduct from this - de Q-map of the unit-matrix  times a constante (scale):\n\n
 * this -= scale* Q(1)
 * @param scale the constant
 */
void TPM::min_qunit(double scale){

   double q = 1.0 + (M - 2*N)*(M - 1.0)/(N*(N - 1.0));

   scale *= q;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= scale;

}

/**
 * calculate the trace of one pair of sp indices of a DPM an put in (*this):\n\n
 * TPM(a,b,d,e) = sum_{c} DPM(a,b,c,d,e,c)
 * @param dpm input DPM
 */
void TPM::bar(DPM &dpm){

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < M;++l)
            (*this)(i,j) += dpm(a,b,l,c,d,l);

      }
   }

   this->symmetrize();

}

/**
 * map a DPM (dpm) on a TPM (*this) with a T1 map, (Q-like map), watch out for the inverse
 * up map, when M = 2*N it is singular! So don't use it!:
 * @param option = +1 T1_down , =-1 inverse T1_up
 * @param dpm The input DPM
 */
void TPM::T(int option,DPM &dpm){

   TPM tpm(M,N);
   tpm.bar(dpm);

   if(option == 1){

      double a = 1;
      double b = 1.0/(3.0*N*(N - 1.0));
      double c = 0.5/(N - 1.0);

      this->Q(1,a,b,c,tpm);

   }
   else{//option == -1

      double a = M - 4.0;
      double b = (M - N - 2.0)/(N*(N - 1.0));
      double c = (M - N - 2.0)/(N - 1.0);

      this->Q(-1,a,b,c,tpm);

   }

}

/**
 * Map a PPHM (pphm) object on a TPM (*this) object by tracing one pair of indices from the pphm (for more info, see primal_dual.pdf)
 * @param pphm input PPHM
 */
void TPM::bar(PPHM &pphm){

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < M;++l)
            (*this)(i,j) += pphm(a,b,l,c,d,l);

      }
   }

   this->symmetrize();

}

/**
 * Map a PPHM (pphm) onto a TPM object (*this) with a T2 down map, see primal_dual.pdf for more information
 * @param pphm input PPHM
 */
void TPM::T(PPHM &pphm){

   //first make some necessary derivate matrices of pphm
   TPM bar(M,N);
   bar.bar(pphm);

   PHM phm(M,N);
   phm.bar(pphm);

   //watch out, scaling for spm is not the usual!
   SPM spm(M,N);

   for(int a = 0;a < M;++a)
      for(int b = a;b < M;++b){

         spm(a,b) = 0;

         for(int c = 0;c < M;++c)
            spm(a,b) += phm(c,a,c,b);

         spm(a,b) *= 0.5/(N - 1.0);

      }

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         //first the tp part:
         (*this)(i,j) = bar(i,j);

         //then the ph part:
         (*this)(i,j) -= phm(d,a,b,c) - phm(d,b,a,c) - phm(c,a,b,d) + phm(c,b,a,d);

         //finaly the three sp parts:
         if(b == d)
            (*this)(i,j) += spm(a,c);

         if(b == c)
            (*this)(i,j) -= spm(a,d);

         if(a == c)
            (*this)(i,j) += spm(b,d);

      }
   }

   this->symmetrize();

}

/**
 * Collaps a SUP matrix S onto a TPM matrix like this:\n\n
 * sum_i Tr (S u^i)f^i = this
 * @param option = 0, project onto full symmetric matrix space, = 1 project onto traceless symmetric matrix space
 * @param S input SUP
 */
void TPM::collaps(int option,SUP &S){

   *this = S.tpm(0);

   TPM hulp(M,N);

   hulp.Q(1,S.tpm(1));

   *this += hulp;

#ifdef __G_CON

   hulp.G(1,S.phm());

   *this += hulp;

#endif

#ifdef __T1_CON

   hulp.T(1,S.dpm());

   *this += hulp;

#endif

#ifdef __T2_CON
   
   hulp.T(S.pphm());

   *this += hulp;

#endif

   if(option == 1)
      this->proj_Tr();

}

/**
 * Print TPM matrix to file called filename
 * @param filename char containing the name and location of the file
 */
void TPM::out(const char *filename){

   this->Matrix::out(filename);

   ofstream output;

   output.precision(10);

   output.open(filename,ios::app);

   output << M << "\t" << N << endl;

}

/** 
 * Construct the pairing hamiltonian with a single particle spectrum
 * @param pair_coupling The strenght of the pairing interaction
 */
void TPM::sp_pairing(double pair_coupling){

   double *E = new double [M/2];

   //single particle spectrum
   for(int a = -M/2;a < 0;++a)
      E[M/2 + a] = (double) a*2/M;

   double *x = new double [M/2];

   //pairing interaction term
   for(int a = 0;a < M/2;++a)
      x[a] = a;//1.0;

   //normeren op 1/2
   double ward = 0.0;

   for(int i = 0;i < M/2;++i)
      ward += x[i]*x[i];

   ward *= 2.0;

   for(int a = 0;a < M/2;++a)
      x[a] /= std::sqrt(ward);

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      (*this)(i,i) = (E[a/2] + E[b/2])/(N - 1.0);

      if(a/2 == b/2)
         (*this)(i,i) -= 2.0*pair_coupling*x[a/2]*x[a/2];

      for(int j = i + 1;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         if(a/2 == b/2 && c/2 == d/2)
            (*this)(i,j) = -2.0*pair_coupling*x[a/2]*x[c/2];
         else
            (*this)(i,j) = 0.0;

      }

   }

   this->symmetrize();

   delete [] E;

}

void TPM::in_sp(const char *filename){

   ifstream input(filename);

   double value;

   int a,b,c,d;

   int i,j;

   while(input >> a >> b >> c >> d >> value){

      i = s2t[a][b];
      j = s2t[c][d];

      (*this)(i,j) = value;

   }

   this->symmetrize();

}

/* vim: set ts=3 sw=3 expandtab :*/
