#ifndef LINCON_H
#define LINCON_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class TPM;
class SPM;

/**
 * @author Brecht Verstichel
 * @date 08-12-2010\n\n
 * This is a class that contains the information about a single linear constraint.
 */

class LinCon{

   /**
    * output stream operator overloaded, will print the constraint matrix, the value of the minimal projection, and the current projection.
    * @param output The stream to which you are writing (e.g. cout)
    * @param lc_p the LinCon object you want to print
    */
   friend ostream &operator<<(ostream &output,const LinCon &lc_p);

   public:

      //constructor
      LinCon(int,int);

      //copy constructor
      LinCon(const LinCon &);

      //destructor
      virtual ~LinCon();

      const TPM &gI() const;

      const SPM &gI_bar() const;

      double gI_tr() const;

      double gi() const;

      void sI(const TPM &);

      void si(double);

      int gM() const;

      int gN() const;

      void diag_T(int);

      void spincon(double);

      void fill_Random();

      void in(const char *filename);

   private:

      //!Traceless Constraint matrix: Watch out, shifted with i_c unity, so that Tr Gamma I_c > 0
      TPM *I_c;

      //!Partially traced, traceless constraint matrix
      SPM *I_c_bar;

      //!scaled trace of the contraint
      double I_c_tr;

      //!minimal projection on the constraint matrix, such that Tr(Gamma I_C) geq i_c
      double i_c;

      int M;

      int N;

};

#endif
