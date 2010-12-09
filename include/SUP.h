#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "TPM.h"
#include "PHM.h"
#include "DPM.h"
#include "PPHM.h"

//definitions:
#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#endif

#ifdef PQGT1

#define __Q_CON
#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __Q_CON
#define __G_CON
#define __T2_CON

#endif

#ifdef PQGT

#define __Q_CON
#define __G_CON
#define __T1_CON
#define __T2_CON

#endif

class EIG;

/**
 * @author Brecht Verstichel
 * @date 09-03-2010\n\n
 * This class, SUP is a blockmatrix over the carrierspace's of active N-representability conditions. 
 * This class contains two TPM objects, and if compiled with the right option a PHM or DPM object, 
 * You have to remember that these matrices are independent of each other (by which I mean that TPM::Q(SUP_PQ::tpm (0))
 * is not neccesarily equal to SUP_PQ::tpm (1)) etc. .
 */
class SUP{
  
   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << sup_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << sup_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param SZ_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,SUP &SZ_p);

   public:

      //constructor
      SUP(int M,int N);

      //copy constructor
      SUP(SUP &);

      //destructor
      ~SUP();

      //overload += operator
      SUP &operator+=(SUP &);

      //overload -= operator
      SUP &operator-=(SUP &);

      //overload equality operator
      SUP &operator=(SUP &);

      //overload equality operator
      SUP &operator=(double );

      TPM &tpm(int i);

      //initialiseer S
      void init_S();

      //initialiseer Z
      void init_Z(double alpha,TPM &ham,SUP &u_0);

      int gN();

      int gM();

      int gn_tp();

      int gdim();

      double ddot(SUP &);

      void invert();

      void dscal(double alpha);

      void proj_U();

      void proj_C(TPM &);

      //maak de matrix D, nodig voor de hessiaan van het stelsel
      void D(SUP &S,SUP &Z);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(SUP &,SUP &);

      void daxpy(double alpha,SUP &);

      double trace();

      double U_trace();

      void proj_C();

      SUP &mprod(SUP &,SUP &);

      void fill(TPM &);

      void fill();

      int solve(SUP &B,SUP &D);

      void H(SUP &B,SUP &D);

      void proj_U_Tr();

      double U_norm();

      double center_dev(SUP &Z);

      double line_search(SUP &DZ,SUP &S,SUP &Z,double max_dev);

      void fill_Random();

      void sep_pm(SUP &p,SUP &m);

#ifdef __G_CON

      PHM &phm();

      int gn_ph();

#endif

#ifdef __T1_CON
      
      DPM &dpm();

      int gn_dp();

#endif

#ifdef __T2_CON

      PPHM &pphm();

      int gn_pph();

#endif

   private:

      //!double pointer of TPM's, will contain the P and Q block of the SUP in the first and second block.
      TPM **SZ_tp;

      //!number of sp orbitals
      int M;

      //!nr of particles
      int N;

      //!dimension of tp space
      int n_tp;

      //!total dimension of the SUP matrix
      int dim;

#ifdef __G_CON

      //!pointer to the particle hole matrix
      PHM *SZ_ph;

      //!dimenson of particle hole space
      int n_ph;

#endif

#ifdef __T1_CON
      
      //!pointer tot he three particles matrix DPM
      DPM *SZ_dp;

      //!dimension of three particle space
      int n_dp;

#endif

#ifdef __T2_CON

      //!pointer tot he three particles matrix DPM
      PPHM *SZ_pph;

      //!dimension of three particle space
      int n_pph;

#endif

};

#endif
