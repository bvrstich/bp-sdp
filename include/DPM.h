#ifndef DPM_H
#define DPM_H

#include <iostream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, DPM, is a class written for three-particle matrices (name comes from drie-particle matrix). It is written specially for the T_1 condition. 
 * It inherits all the functions from its mother class Matrix, some special member functions and two lists that give the relationship between the dp (three-particle) and the sp basis.
 */
class DPM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << dpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << dpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param dpm_p the DPM you want to print
    */
   friend ostream &operator<<(ostream &output,DPM &dpm_p);

   public:
      
      //constructor
      DPM(int M,int N);

      //copy constructor
      DPM(const DPM &);

      //destructor
      virtual ~DPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int f) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //geef dim terug
      int gn() const;

      //generalized T1 map
      void T(double,double,double,const TPM &);

      //maak een DPM van een TPM via de T1 conditie
      void T(int option,const TPM &);

      //maak een DPM van een TPM via de hat functie
      void hat(const TPM &);

      //deduct scale times T1(1) matrix
      void min_tunit(double scale);

      //input from file with sp indices
      void in_sp(const char *);

   private:

      //!static counter that counts the number of DPM objects running in the program
      static int counter;

      //!static list of dimension [n_dp][3] that takes in a dp index i and returns three sp indices: a = dp2s[i][0], b = dp2s[i][1] and c = dp2s[i][2]
      static int **dp2s;

      //!static list of dimension [M][M][M] that takes three sp indices a,b and c and returns a dp index i: i = s2dp[a][b][c]
      static int ***s2dp;

      //!nr of particles
      int N;

      //!dimension of sp hilbert space
      int M;

      //!dimension of dp space
      int n;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
