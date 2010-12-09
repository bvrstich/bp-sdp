#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * This class SPM was written for single particle matrices. It inherits from the class Matrix and expands it with
 * specific memberfunction and a knowledge of the nr of sp orbitals and particles.
 */

class SPM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << spm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << spm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de SPM you want to print
    */
   friend ostream &operator<<(ostream &output,SPM &spm_p);

   public:
      
      //constructor
      SPM(int M,int N);

      //copy constructor
      SPM(SPM &);

      //destructor
      virtual ~SPM();

      using Matrix::operator=;

      int gN();

      int gM();

      /**
       * constructs a SPM from a TPM or a PHM. Definition of these functions are in different notes.
       * Actually this function is defined as bar * scale.
       * @param scale The factor that multiplies the bar(MT), e.g. 1/(N - 1) for a normal single particle density matrix
       * @param MT PHM or TPM inputmatrix.
       */
      template<class MatrixType>
         void constr(double scale,MatrixType &MT){

            for(int a = 0;a < M;++a)
               for(int b = a;b < M;++b){

                  (*this)(a,b) = 0.0;

                  for(int l = 0;l < M;++l)
                     (*this)(a,b) += MT(a,l,b,l);

                  (*this)(a,b) *= scale;

               }

            this->symmetrize();

         }

      /**
       * Constructor of an SPM object with initialisation by means of the function SPM::constr
       * @param scale The factor that multiplies the bar(MT), e.g. 1/(N - 1) for a normal single particle density matrix
       * @param MT the PHM or TPM inputmatrix.
       */
      template<class MatrixType>
         SPM(double scale,MatrixType &MT) : Matrix(MT.gM()) {

            this->M = MT.gM();
            this->N = MT.gN();

            this->constr(scale,MT);

         }

      /**
       * construeert een SPM uit een TPM of PHM. \n\n
       * SPM(a,c) = sum_b MT(a,b,c,b)\n\n
       * @param MT the PHM or TPM inputmatrix.
       */
      template<class MatrixType>
         void bar(MatrixType &MT){

            for(int a = 0;a < M;++a)
               for(int b = a;b < M;++b){

                  (*this)(a,b) = 0.0;

                  for(int l = 0;l < M;++l)
                     (*this)(a,b) += MT(a,l,b,l);

               }

            this->symmetrize();

         }

   private:

      //!dimension of single particle space
      int M;

      //!nr of particles
      int N;

};

#endif
