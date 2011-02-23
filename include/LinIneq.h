#ifndef LININEQ_H
#define LININEQ_H

#include <iostream>
#include <cstdlib>

using std::ostream;

#include "LinCon.h"

/**
 * @author Brecht Verstichel
 * @date 08-12-2010\n\n
 * This is a class that contains the information about all the constraints active at a certain point in the program.
 */

class LinIneq{

   /**
    * output stream operator overloaded, will print all of the LinCon objects
    * @param output The stream to which you are writing (e.g. cout)
    * @param li_p the LinIneq object you want to print
    */
   friend ostream &operator<<(ostream &output,const LinIneq &li_p);

   public:

      //constructor
      LinIneq(int,int);

      //copy constructor
      LinIneq(const LinIneq &);

      //destructor
      virtual ~LinIneq();

      LinIneq &operator=(const LinIneq &);

      LinIneq &operator+=(const LinIneq &);

      LinIneq &operator-=(const LinIneq &);

      LinIneq &operator=(double);

      int gnr() const;

      //easy to access the LinCon objects
      const LinCon &operator[](int i) const;

      //easy to access and change the LinCon objects
      LinCon &operator[](int i);

      int gM() const;

      int gN() const;

      void fill(const TPM &);

      double gproj(int) const;

      double *gproj();

      double gproj_bar(int) const;

      double *gproj_bar();

      double gtr() const;

      static void init(int,int,int);

      static void clean();

      static void constr_overlap(int,int);

      double ga() const;

      double gc() const;

      double alpha(int) const;

      void fill_Random();

      double ddot(const LinIneq &) const;

      void invert();

      void dcal(double);

      void sqrt(int);

      void L_map(const LinIneq &,const LinIneq &);

      LinIneq &daxpy(double,const LinIneq &);

      void dscal(double);

      static void print_coef();

   private:

      //!LinCon array containing the different LinCon objects
      static LinCon **li;

      //!counter of the nr of objects in the program
      static int counter;

      //!nr of linear constraints
      static int nr;

      //!variables of the overlapmatrix-map
      static double a,c;

      //!coefficient matrix of the overlap matrix map system of linear equations
      static double *coef;

      //!array containing the projections on the constraint
      double *proj;

      //!array containing the barred projections onto the barred constraints
      double *proj_bar;

      //!the trace of the input TPM object.
      double tr;

      //!nr of sp orbitals
      int M;

      //!nr of particles
      int N;

};

#endif
