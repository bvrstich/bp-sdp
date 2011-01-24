#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

class SPM;
class SUP;
class PHM;
class DPM;
class PPHM;

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class TPM is a class written for two particle matrices, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */

class TPM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,TPM &tpm_p);

   public:
      
      //constructor
      TPM(int M,int N);

      //copy constructor
      TPM(TPM &);

      //file constructor
      TPM(const char *);

      //destructor
      virtual ~TPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d) const;

      //geef N terug
      int gN();

      //geef N terug
      int gM();

      //geef n terug
      int gn();

      void hubbard(int option,double U);

      //Q afbeelding en zijn inverse
      void Q(int option,TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,TPM &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int option,TPM &);

      void unit();

      void proj_Tr();

      //de hessiaan afbeelding:
      void H(TPM &b,SUP &D);

      //los het stelsel op
      int solve(TPM &b,SUP &D);

      //de G down en inverse G up
      void G(int option,PHM &);

      //trace one pair of indices of DPM
      void bar(DPM &);

      //trace one pair of indices of PPHM
      void bar(PPHM &);

      //T1 down
      void T(int option,DPM &);

      //T2 down
      void T(PPHM &);

      void min_unit(double scale);

      void min_qunit(double scale);

      void collaps(int option,SUP &);

      void out(const char *filename);

      void sp_pairing(double );

      void in_sp(const char *);

   private:

      //!static list of dimension [n_tp][2] that takes in a tp index i and returns two sp indices: a = t2s[i][0] and b = t2s[i][1]
      static int **t2s;

      //!static list of dimension [M][M] that takes two sp indices a,b and returns a tp index i: i = s2t[a][b]
      static int **s2t;

      //!static counter that counts the number of TPM objects running in the program
      static int counter;

      //!nr of particles
      int N;

      //!dimension of sp hilbert space
      int M;

      //!dimension of tp hilbert space and of the matrix
      int n;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
