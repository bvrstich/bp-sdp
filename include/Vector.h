#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;

#include "lapack.h"

template<class MatrixType>
class Vector;

/**
 * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
 * ifstream object and type:\n\n
 * object << vector << endl;\n\n
 * For output onto the screen type: \n\n
 * cout << vector << endl;\n\n
 * @param output The stream to which you are writing (e.g. cout)
 * @param vector_p de Vector you want to print
 */
template<class MatrixType>
ostream &operator<<(ostream &output,Vector<MatrixType> &vector_p);

/**
 * @author Brecht Verstichel
 * @date 15-04-2010\n\n
 * This is a class written for vectors. It will contain the eigenvalues of the TPM, etc. Matrices. It is a template class,
 * corresponding to the different VectorType's that can be put in, it will automatically get the right dimension.
 * It is a wrapper around a pointer and redefines much used lapack and blas routines as memberfunctions
 */
template<class MatrixType>
class Vector{

   public:

      //construct with as input a MatrixType
      Vector(MatrixType& );

      //copy constructor
      Vector(Vector<MatrixType> &);

      //destructor
      virtual ~Vector();

      //overload equality operator
      Vector &operator=(Vector<MatrixType> &);

      Vector &operator=(double );

      //overload += operator
      Vector &operator+=(Vector<MatrixType> &);

      //overload -= operator
      Vector &operator-=(Vector<MatrixType> &);

      Vector &daxpy(double alpha,Vector<MatrixType> &);

      Vector &operator/=(double );

      //easy to change the numbers
      double &operator[](int i);

      //easy to access the numbers
      double operator[](int i) const;

      //get the pointer to the vector
      double *gVector();

      int gn();

      double sum();

      double log_product();

      void diagonalize(MatrixType &);

      double ddot(Vector<MatrixType> &);

      void dscal(double alpha);

      double min();
      
      double max();

      double centerpot(double );

   private:

      //!pointer of doubles, contains the numbers, the vector
      double *vector;

      //!dimension of the vector
      int n;

};

/**
 * Construct and initialize the Vector object by diagonalizing a Matrix object
 */
template<class MatrixType>
Vector<MatrixType>::Vector(MatrixType &MT){

   //allocate
   this->n = MT.gn();

   vector = new double [n];

   //initialize
   char jobz = 'V';
   char uplo = 'U';

   int lwork = 3*n - 1;

   double *work = new double [lwork];

   int info;

   dsyev_(&jobz,&uplo,&n,(MT.gMatrix())[0],&n,vector,work,&lwork,&info);

   delete [] work;

}

/**
 * copy constructor 
 * @param vec_copy The vector you want to be copied into the object you are constructing, make sure that it is an allocated and filled vector!
 */
template<class MatrixType>
Vector<MatrixType>::Vector(Vector<MatrixType> &vec_copy){

   this->n = vec_copy.n;

   vector = new double [n];

   int inc = 1;

   dcopy_(&n,vec_copy.vector,&inc,vector,&inc);

}

/**
 * Destructor
 */
template<class MatrixType>
Vector<MatrixType>::~Vector(){

   delete [] vector;

}

/**
 * overload the equality operator
 * @param vector_copy The vector you want to be copied into this
 */
template<class MatrixType>
Vector<MatrixType> &Vector<MatrixType>::operator=(Vector<MatrixType> &vector_copy){

   int inc = 1;

   dcopy_(&n,vector_copy.vector,&inc,vector,&inc);

   return *this;

}

/**
 * Make all the number in your vector equal to the number a, e.g. usefull for initialization (Vector M = 0)
 * @param a the number
 */
template<class MatrixType>
Vector<MatrixType> &Vector<MatrixType>::operator=(double a){

   for(int i = 0;i < n;++i)
      vector[i] = a;

   return *this;

}

/**
 * overload the += operator for matrices
 * @param vector_pl The vector you want to add to this
 */
template<class MatrixType>
Vector<MatrixType> &Vector<MatrixType>::operator+=(Vector<MatrixType> &vector_pl){

   int inc = 1;
   double alpha = 1.0;

   daxpy_(&n,&alpha,vector_pl.vector,&inc,vector,&inc);

   return *this;

}

/**
 * overload the -= operator for matrices
 * @param vector_pl The vector you want to deduct from this
 */
template<class MatrixType>
Vector<MatrixType> &Vector<MatrixType>::operator-=(Vector<MatrixType> &vector_pl){

   int inc = 1;
   double alpha = -1.0;

   daxpy_(&n,&alpha,vector_pl.vector,&inc,vector,&inc);

   return *this;

}

/**
 * add the vector vector_pl times the constant alpha to this
 * @param alpha the constant to multiply the vector_pl with
 * @param vector_pl the Vector to be multiplied by alpha and added to this
 */
template<class MatrixType>
Vector<MatrixType> &Vector<MatrixType>::daxpy(double alpha,Vector<MatrixType> &vector_pl){

   int inc = 1;

   daxpy_(&n,&alpha,vector_pl.vector,&inc,vector,&inc);

   return *this;

}

/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your vector through
 */
template<class MatrixType>
Vector<MatrixType> &Vector<MatrixType>::operator/=(double c){

   int inc = 1;

   double alpha = 1.0/c;

   dscal_(&n,&alpha,vector,&inc);

   return *this;

}

/**
 * write access to your vector, change the number on index i
 * @param i row number
 * @return the entry on place i
 */
template<class MatrixType>
double &Vector<MatrixType>::operator[](int i){

   return vector[i];

}

/**
 * read access to your vector, view the number on index i
 * @param i row number
 * @return the entry on place i
 */
template<class MatrixType>
double Vector<MatrixType>::operator[](int i) const {

   return vector[i];

}

/**
 * @return the underlying pointer to vector, useful for mkl and lapack applications
 */
template<class MatrixType>
double *Vector<MatrixType>::gVector(){

   return vector;

}

/**
 * @return the dimension of the vector
 */
template<class MatrixType>
int Vector<MatrixType>::gn(){

   return n;

}

/**
 * @return the sum of all the elements in the vector
 */
template<class MatrixType>
double Vector<MatrixType>::sum(){

   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += vector[i];

   return ward;

}

/**
 * @return the logarithm of the product of all the elements in the vector (so the sum of all the logarithms)
 */
template<class MatrixType>
double Vector<MatrixType>::log_product(){

   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += log(vector[i]);

   return ward;

}

/**
 * diagonalize a MatrixType object when the memory has allready been allocated, make sure the MT object has the same dimension as the Vector.
 * @param MT the MatrixType object that will be diagonalized
 */
template<class MatrixType>
void Vector<MatrixType>::diagonalize(MatrixType &MT){

   //initialize
   char jobz = 'V';
   char uplo = 'U';

   int lwork = 3*n - 1;

   double *work = new double [lwork];

   int info;

   dsyev_(&jobz,&uplo,&n,(MT.gMatrix())[0],&n,vector,work,&lwork,&info);

   delete [] work;

}

/**
 * @return inproduct of (*this) vector with vector_i
 * @param vector_i input vector
 */
template<class MatrixType>
double Vector<MatrixType>::ddot(Vector &vector_i){

   int inc = 1;

   return ddot_(&n,vector,&inc,vector_i.vector,&inc);

}

/**
 * Scale the vector (*this) with parameter alpha
 * @param alpha scalefactor
 */
template<class MatrixType>
void Vector<MatrixType>::dscal(double alpha){

   int inc = 1;

   dscal_(&n,&alpha,vector,&inc);

}

template<class MatrixType>
ostream &operator<<(ostream &output,Vector<MatrixType> &vector_p){

   for(int i = 0;i < vector_p.gn();++i)
      output << i << "\t" << vector_p[i] << std::endl;

   return output;

}

/**
 * @return the minimal element present in this Vector object.
 * watch out, only works when Vector is filled with the eigenvalues of a diagonalized Matrix object
 */
template<class MatrixType>
double Vector<MatrixType>::min(){

   return vector[0];

}

/**
 * @return the maximal element present in this Vector object.
 * watch out, only works when Vector is filled with the eigenvalues of a diagonalized Matrix object
 */
template<class MatrixType>
double Vector<MatrixType>::max(){

   return vector[n - 1];

}

template<class MatrixType>
double Vector<MatrixType>::centerpot(double alpha){

   double ward = 0.0;

   for(int i = 0;i < n;++i)
      ward += log(1.0 + alpha*vector[i]);

   return ward;

}

#endif
