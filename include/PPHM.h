#ifndef PPHM_H
#define PPHM_H

#include <iostream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 18-03-2010\n\n
 * This class, PPHM, is a class written for two-particle-one-hole matrices. It is written specially for the T_2 condition. 
 * It inherits all the functions from its mother class Matrix, some special member functions and two lists that give the relationship between the pph (two-particle one hole) and the sp basis.
 */
class PPHM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << pphm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << pphm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param pphm_p the PPHM you want to print
    */
   friend ostream &operator<<(ostream &output,PPHM &pphm_p);

   public:
      
      //constructor
      PPHM(int M,int N);

      //copy constructor
      PPHM(PPHM &);

      //destructor
      virtual ~PPHM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int f) const;

      //geef N terug
      int gN();

      //geef M terug
      int gM();

      //geef dim terug
      int gn();

      //maak een PPHM van een TPM via de T2 conditie
      void T(int option,TPM &);

      void min_tunit(double );

      double skew_trace();

   private:

      //!static counter that counts the number of PPHM objects running in the program
      static int counter;

      //!static list of dimension [n_pph][3] that takes in a pph index i and returns three sp indices: a = pph2s[i][0], b = pph2s[i][1] and c = pph2s[i][2]
      static int **pph2s;

      //!static list of dimension [M][M][M] that takes three sp indices a,b and c and returns a pph index i: i = s2pph[a][b][c]
      static int ***s2pph;

      //!nr of particles
      int N;

      //!dimension of sp hilbert space
      int M;

      //!dimension of pph space
      int n;

};

#endif
