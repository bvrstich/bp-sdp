#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <cstdlib>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 28-02-2011\n\n
 * This is a class written for the 2 dimension hubbard model, it translates the physical degrees of freedom of the
 * 2D hubbard model: (x,y,sigma) to the one particle basis of the regular program.
 */

class Hamiltonian{

   public:

      //!initializes the lists.
      static void init(int);

      //!clears the lists;
      static void clear();

      //!print the list
      static void print();

      static void construct_T(Matrix &);

   private:
      
      //!the dimension of the lattice
      static int L;

      //!the dimension of the sp-hilbert space
      static int M;

      //!static list that translates the three indices of the 2D hubbard model to one sp index.
      static int ***xys_alpha;

      //!static list that translates the sp index alpha to the three physical indices of the 2D hubbard model.
      static int **alpha_xys;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
