#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::cout;
using std::endl;

#include "include.h"

/**
 * standard constructor\n
 * Allocates two TPM matrices and optionally a PHM, DPM or PPHM matrix. And of course also the LinIneq object.
 * @param M number of sp orbitals
 * @param N number of particles
 */
SUP::SUP(int M,int N){

   this->M = M;
   this->N = N;
   this->n_tp = M*(M - 1)/2;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   this->dim = 2*n_tp;

#ifdef __G_CON

   this->n_ph = M*M;
   
   SZ_ph = new PHM(M,N);

   dim += n_ph;

#endif

#ifdef __T1_CON
   
   this->n_dp = M*(M - 1)*(M - 2)/6;

   SZ_dp = new DPM(M,N);

   dim += n_dp;

#endif

#ifdef __T2_CON

   this->n_pph = M*M*(M - 1)/2;

   SZ_pph = new PPHM(M,N);

   dim += n_pph;

#endif

   li = new LinIneq(M,N);

}

/**
 * standard constructor\n
 * Allocates two TPM matrices and optionally a PHM, DPM or PPHM matrix, and a LinIneq object, then copies the content of
 * input SUP SZ_c into it.
 * @param SZ_c input SUP
 */
SUP::SUP(const SUP &SZ_c)
{
   this->M = SZ_c.M;
   this->N = SZ_c.N;
   this->n_tp = SZ_c.n_tp;
   this->dim = 2*n_tp;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

#ifdef __G_CON

   this->n_ph = M*M;

   dim += n_ph;
   
   SZ_ph = new PHM(M,N);

   *SZ_ph = *SZ_c.SZ_ph;

#endif

#ifdef __T1_CON
   
   this->n_dp = M*(M - 1)*(M - 2)/6;

   SZ_dp = new DPM(M,N);

   dim += n_dp;

   *SZ_dp = *SZ_c.SZ_dp;

#endif

#ifdef __T2_CON

   this->n_pph = M*M*(M - 1)/2;

   SZ_pph = new PPHM(M,N);

   dim += n_pph;

   *SZ_pph = *SZ_c.SZ_pph;

#endif

   li = new LinIneq(SZ_c.gli());

   dim += li->gnr();

}

/**
 * Destructor
 */
SUP::~SUP(){

   for(int i = 0;i < 2;++i)
      delete SZ_tp[i];

   delete [] SZ_tp;

#ifdef __G_CON
   
   delete SZ_ph;

#endif

#ifdef __T1_CON

   delete SZ_dp;

#endif

#ifdef __T2_CON
   
   delete SZ_pph;

#endif

   delete li;

}

/**
 * Overload += operator
 * @param SZ_pl The SUP matrix that has to be added to this
 */
SUP &SUP::operator+=(const SUP &SZ_pl)
{
   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) += (*SZ_pl.SZ_tp[i]);

#ifdef __G_CON
   
   (*SZ_ph) += (*SZ_pl.SZ_ph);

#endif

#ifdef __T1_CON

   (*SZ_dp) += (*SZ_pl.SZ_dp);

#endif

#ifdef __T2_CON

   (*SZ_pph) += (*SZ_pl.SZ_pph);

#endif

   (*li) += SZ_pl.gli();

   return *this;

}

/**
 * Overload -= operator
 * @param SZ_pl The SUP that will be deducted from this
 */
SUP &SUP::operator-=(const SUP &SZ_pl)
{
   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) -= (*SZ_pl.SZ_tp[i]);

#ifdef __G_CON
   
   (*SZ_ph) -= (*SZ_pl.SZ_ph);

#endif

#ifdef __T1_CON

   (*SZ_dp) -= (*SZ_pl.SZ_dp);

#endif

#ifdef __T2_CON

   (*SZ_pph) -= (*SZ_pl.SZ_pph);

#endif

   (*li) -= SZ_pl.gli();

   return *this;

}

/**
 * Overload equality operator, copy SZ_c into this
 * @param SZ_c SUP_PQ to be copied into this
 */
SUP &SUP::operator=(const SUP &SZ_c)
{
   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

#ifdef __G_CON

   (*SZ_ph) = (*SZ_c.SZ_ph);

#endif

#ifdef __T1_CON

   (*SZ_dp) = (*SZ_c.SZ_dp);

#endif

#ifdef __T2_CON

   (*SZ_pph) = (*SZ_c.SZ_pph);

#endif

   (*li) = SZ_c.gli();

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. SZ = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP &SUP::operator=(const double &a)
{
   (*SZ_tp[0]) = a;
   (*SZ_tp[1]) = a;

#ifdef __G_CON

   (*SZ_ph) = a;

#endif

#ifdef __T1_CON

   (*SZ_dp) = a;

#endif

#ifdef __T2_CON

   (*SZ_pph) = a;

#endif

   (*li) = a;

   return *this;

}

/**
 * The const version
 * @param i which block you want to have the pointer to.
 * @return pointer to the individual TPM blocks: SZ_tp[i]
 */
const TPM &SUP::tpm(int i) const
{
   return *SZ_tp[i];
}

/**
 * @param i which block you want to have the pointer to.
 * @return pointer to the individual TPM blocks: SZ_tp[i]
 */
TPM &SUP::tpm(int i)
{
   return *SZ_tp[i];
}

#ifdef __G_CON

/**
 * The const version
 * @return pointer to the PHM block: SZ_ph
 */
const PHM &SUP::phm() const
{
   return *SZ_ph;
}

/**
 * @return pointer to the PHM block: SZ_ph
 */
PHM &SUP::phm()
{
   return *SZ_ph;
}

#endif

#ifdef __T1_CON

/**
 * The const version
 * @return pointer to the DPM block: SZ_dp
 */
const DPM &SUP::dpm() const
{
   return *SZ_dp;
}

/**
 * @return pointer to the DPM block: SZ_dp
 */
DPM &SUP::dpm()
{
   return *SZ_dp;
}

#endif

#ifdef __T2_CON

/**
 * The const version
 * @return pointer to the PPHM block: SZ_pph
 */
const PPHM &SUP::pphm() const
{
   return *SZ_pph;
}

/**
 * @return pointer to the PPHM block: SZ_pph
 */
PPHM &SUP::pphm()
{
   return *SZ_pph;
}

#endif

/**
 * The const version
 * @return pointer to the LinIneq object li.
 */
const LinIneq &SUP::gli() const{

   return *li;

}

/**
 * @return pointer to the LinIneq object li.
 */
LinIneq &SUP::gli(){

   return *li;

}

/**
 * Initialization of the SUP matrix S, is just u^0: see primal_dual.pdf for more information
 */
void SUP::init_S(){

   (*SZ_tp[0]).unit();

   this->fill();

}

ostream &operator<<(ostream &output,const SUP &SZ_p){

   output << (*SZ_p.SZ_tp[0]) << std::endl;
   output << (*SZ_p.SZ_tp[1]);

#ifdef __G_CON

   output << std::endl;
   output << (*SZ_p.SZ_ph);

#endif

#ifdef __T1_CON

   output << std::endl;
   output << (*SZ_p.SZ_dp);

#endif

#ifdef __T2_CON

   output << std::endl;
   output << (*SZ_p.SZ_pph);

#endif

   output << std::endl;
   output << SZ_p.gli();

   return output;

}

/**
 * Fill the SUP matrix with random elements, watch out, not necessarily positive definite
 */
void SUP::fill_Random(){

   SZ_tp[0]->fill_Random();
   SZ_tp[1]->fill_Random();

#ifdef __G_CON

   SZ_ph->fill_Random();

#endif

#ifdef __T1_CON

   SZ_dp->fill_Random();

#endif

#ifdef __T2_CON

   SZ_pph->fill_Random();

#endif

   li->fill_Random();

}

/**
 * Initialisation for dual SUP matrix Z, see primal_dual.pdf for info.
 */
void SUP::init_Z(double alpha,const TPM &ham,const SUP &u_0)
{
   this->fill_Random();

   //nog een eenheidsmatrix maal constante bijtellen zodat Z sterk positief definiet is:
   this->daxpy(alpha,u_0); 

   //en projecteren!
   this->proj_C(ham);

}

/**
 * @return number of particles
 */
int SUP::gN() const
{
   return N;
}

/**
 * @return dimension of sp space
 */
int SUP::gM() const
{
   return M;
}

/**
 * @return dimension of tp space
 */
int SUP::gn_tp() const
{
   return n_tp;
}

#ifdef __G_CON

/**
 * @return dimension of ph space
 */
int SUP::gn_ph() const
{
   return n_ph;
}

#endif

#ifdef __T1_CON

/**
 * @return dimension of dp space
 */
int SUP::gn_dp() const
{
   return n_dp;
}

#endif

#ifdef __T2_CON

/**
 * @return dimension of pph space
 */
int SUP::gn_pph() const
{
   return n_pph;
}

#endif

/**
 * @return total dimension of SUP (carrier) space
 */
int SUP::gdim() const
{
   return dim;
}

/**
 * @return nr of constraints
 */
int SUP::gnr() const {

   return li->gnr();

}

/**
 * @param SZ_i input SUP_PQ SZ_i
 * @return inproduct between this and input matrix SZ_i, defined as Tr(this SZ_i)
 */
double SUP::ddot(const SUP &SZ_i) const
{
   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->ddot(*SZ_i.SZ_tp[i]);

#ifdef __G_CON
   
   ward += SZ_ph->ddot(*SZ_i.SZ_ph);

#endif

#ifdef __T1_CON

   ward += SZ_dp->ddot(*SZ_i.SZ_dp);

#endif

#ifdef __T2_CON

   ward += SZ_pph->ddot(*SZ_i.SZ_pph);

#endif

   ward += li->ddot(SZ_i.gli());

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 * Makes use of cholesky decomposition, so only positive definite matrices can be used as input!
 */
void SUP::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

#ifdef __G_CON
   
   SZ_ph->invert();

#endif

#ifdef __T1_CON
   
   SZ_dp->invert();

#endif

#ifdef __T2_CON
   
   SZ_pph->invert();

#endif

   li->invert();

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

#ifdef __G_CON
   
   SZ_ph->dscal(alpha);

#endif

#ifdef __T1_CON
   
   SZ_dp->dscal(alpha);

#endif

#ifdef __T2_CON
   
   SZ_pph->dscal(alpha);

#endif

   li->dscal(alpha);

}

/**
 * Orthogonal projection of a general SUP matrix diag[ M M_Q ( M_G M_T1 M_T2 ) ] onto U space: diag[ M_u Q(M_u) ( G(M_u) T1(M_u) T2(M_u) ) ]
 * for more information see primal_dual.pdf
 */
void SUP::proj_U(){
  
   //eerst M_Gamma + Q(M_Q) + ( G(M_G) + T1(M_T1) + T2(M_T2) ) in O stoppen
   TPM O(M,N);

   //collaps onto traceless TPM space
   O.collaps(1,*this);

   //dan de inverse overlapmatrix hierop laten inwerken en in this[0] stoppen
   SZ_tp[0]->S_L(-1,O);

   //fill up the rest with the right maps
   this->fill();

}

/**
 * Project the general SUP matrix (*this) orthogonally onto the linear space for which\n\n
 * Tr(Z u^i) = h^i      with h^i = Tr(tpm f^i)\n\n
 * is valid.
 * @param tpm input TPM (mostly the hamiltonian of the problem)
 */
void SUP::proj_C(const TPM &tpm)
{
   TPM hulp(M,N);

   hulp.collaps(0,*this);

   hulp -= tpm;

   hulp.proj_Tr();

   //Z_res is the orthogonal piece of this that will be deducted,
   //so the piece of this in the U-space - ham
   SUP Z_res(M,N);

   //apply iverse S to it and put it in Z_res.tpm(0)
   (Z_res.tpm(0)).S_L(-1,hulp);

   //and fill it up Johnny
   Z_res.fill();

   *this -= Z_res;

}

/**
 * Construct the D matrix and put it in this, D is the matrix matrix of the hessian, see primal_dual.pdf for more information
 * @param S The primal SUP matrix S
 * @param Z The dual SUP matrix Z
 */
void SUP::D(const SUP &S,const SUP &Z)
{
   //positieve vierkantswortel uit Z
   SUP Z_copy(Z);

   Z_copy.sqrt(1);

   //links en rechts vermenigvuldigen met wortel Z
   SUP hulp(M,N);

   hulp.L_map(Z_copy,S);

   hulp.sqrt(1);

   //negatieve vierkantswortel uit Z
   Z_copy = Z;

   Z_copy.sqrt(-1);

   //en links en rechts hulp vermenigvuldigen met die wortel, en in this steken:
   this->L_map(Z_copy,hulp);

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

#ifdef __G_CON

   SZ_ph->sqrt(option);

#endif

#ifdef __T1_CON

   SZ_dp->sqrt(option);

#endif

#ifdef __T2_CON

   SZ_pph->sqrt(option);

#endif

   li->sqrt(option);

}

/**
 * Multiply symmetric SUP blockmatrix object left en right with symmetric SUP blockmatrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map SUP that will be multiplied to the left en to the right of matrix object
 * @param object central SUP
 */
void SUP::L_map(const SUP &map,const SUP &object)
{
   for(int i = 0;i < 2;++i)
      SZ_tp[i]->L_map(map.tpm(i),object.tpm(i));

#ifdef __G_CON

   SZ_ph->L_map(map.phm(),object.phm());

#endif

#ifdef __T1_CON

   SZ_dp->L_map(map.dpm(),object.dpm());

#endif

#ifdef __T2_CON

   SZ_pph->L_map(map.pphm(),object.pphm());

#endif

   li->L_map(map.gli(),object.gli());

}

/**
 * add the SUP SZ_p times the constant alpha to this
 * @param alpha the constant to multiply the SZ_p with
 * @param SZ_p the SUP to be multiplied by alpha and added to (*this)
 */
void SUP::daxpy(double alpha,const SUP &SZ_p)
{
   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,SZ_p.tpm(i));

#ifdef __G_CON
   
   SZ_ph->daxpy(alpha,SZ_p.phm());

#endif

#ifdef __T1_CON
   
   SZ_dp->daxpy(alpha,SZ_p.dpm());

#endif

#ifdef __T2_CON
   
   SZ_pph->daxpy(alpha,SZ_p.pphm());

#endif

   li->daxpy(alpha,SZ_p.gli());

}

/**
 * Orthogonal projection of a general SUP matrix [ M M_Q ( M_G M_T1 M_T2 ) ] onto the orthogonal complement of the U space (C space)
 * See primal_dual.pdf for more information
 */
void SUP::proj_C(){

   SUP Z_copy(*this);

   //projecteer op de U ruimte
   Z_copy.proj_U();

   //en het orthogonaal complement nemen:
   *this -= Z_copy;

}

/**
 * Fill the SUP matrix (*this) with a TPM matrix like: this = diag[tpm  Q(tpm)  ( G(tpm) T1(tpm) T2(tpm) ) ]
 * @param tpm input TPM
 */
void SUP::fill(const TPM &tpm)
{
   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(1,tpm);

#ifdef __G_CON
   
   SZ_ph->G(1,tpm);

#endif

#ifdef __T1_CON
   
   SZ_dp->T(1,tpm);

#endif

#ifdef __T2_CON
   
   SZ_pph->T(0,tpm);

#endif

   li->fill(tpm);

}

/**
 * fill the SUP matrix with the TPM matrix stored in the first block:\n\n
 * this = diag[this->tpm(0) Q(this->tpm(0)) ( G(this->tpm(0)) T1(this->tpm(0)) T2(this-tpm(0)) ) ]
 */
void SUP::fill(){

   SZ_tp[1]->Q(1,*SZ_tp[0]);

#ifdef __G_CON

   SZ_ph->G(1,*SZ_tp[0]);

#endif 

#ifdef __T1_CON

   SZ_dp->T(1,*SZ_tp[0]);

#endif 

#ifdef __T2_CON

   SZ_pph->T(0,*SZ_tp[0]);

#endif 

   li->fill(*SZ_tp[0]);

}

/**
 * Implementation of the linear conjugate gradient algorithm for the solution of the dual Newton system\n\n
 * H(*this) = B in which H is the dual hessian map
 * @param B right hand side of the equation
 * @param D SUP matrix that defines the structure of the hessian map (the metric) (inverse of the primal Newton equation hessian)
 * @return return the number of iteration required to converge
 */
int SUP::solve(SUP &B,const SUP &D)
{
   SUP HB(M,N);
   HB.H(*this,D);

   B -= HB;

   //de r initialiseren op B - H DZ
   SUP r(B);

   double rr = r.ddot(r);
   double rr_old,ward;

   int cg_iter = 0;

   while(rr > 1.0e-5){

      ++cg_iter;

      HB.H(B,D);

      ward = rr/B.ddot(HB);

      //delta += ward*b
      this->daxpy(ward,B);

      //r -= ward*HB
      r.daxpy(-ward,HB);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      B.dscal(rr/rr_old);

      //dan r er bijtellen
      B += r;

   }
   
   return cg_iter;

}

/**
 * The dual hessian map:\n\n
 * HB = DBD (dus SUP::L_map), projected onto C-space (SUP::proj_C)
 * @param B SUP matrix onto which the hessian works.
 * @param D SUP matrix that defines the structure of the map (metric)
 */
void SUP::H(const SUP &B,const SUP &D)
{
   this->L_map(D,B);

   this->proj_C();

}

/**
 * @return Deviation from the central path measured through the logarithmic potential, it's a measure for
 * the deviation of the product of the primal with the dual matrix (SZ) from the unit matrix.\n
 * Usage of the function: S.center_dev(Z) gives returns the deviation.\n\n
 * (*this) = S = primal matrix of the problem
 * @param Z = dual matrix of the problem
 */
double SUP::center_dev(const SUP &Z) const
{
   SUP sqrt_S(*this);

   sqrt_S.sqrt(1);

   SUP SZ(M,N);
   SZ.L_map(sqrt_S,Z);

   EIG eig(SZ);

   return eig.center_dev();

}

/**
 * Line search function that checks how large a step you can take in a given primal dual predictor direction (DS,DZ), starting from 
 * the current primal dual point (S,Z), before deviating beyond max_dev from the central path.\n\n
 * (*this) = DS --> primal search direction
 * @param DZ dual search direction
 * @param S Current primal point
 * @param Z Current dual point
 * @param max_dev number (double) input by which you can tell the function how far you want to deviate from the central path after the step.
 */
double SUP::line_search(const SUP &DZ,const SUP &S,const SUP &Z,double max_dev) const
{
   //eerst de huidige deviatie van het centraal pad nemen:
   double center_dev = S.center_dev(Z);

   //eigenwaarden zoeken van S^{-1/2} DS S^{-1/2} en Z^{-1/2} DZ Z^{-1/2}

   //kopieer S in de zogeheten wortel:
   SUP wortel(S);

   //maak negatieve vierkantswortel uit S
   wortel.sqrt(-1);

   //de L_map
   SUP hulp(M,N);
   hulp.L_map(wortel,*this);

   //eigenwaarden in eigen_S stoppen
   EIG eigen_S(hulp);

   //nu idem voor Z
   wortel = Z;

   wortel.sqrt(-1);

   hulp.L_map(wortel,DZ);

   EIG eigen_Z(hulp);

   //nog c_S en c_Z uitrekenen:
   double pd_gap = S.ddot(Z);

   //c_S = Tr (DS Z)/Tr (SZ)
   double c_S = this->ddot(Z)/pd_gap;

   //c_Z = Tr (S DZ)/Tr (SZ)
   double c_Z = S.ddot(DZ)/pd_gap;

   //waar zitten de singulariteiten: tot waar mag ik zoeken?
   double a_max = -1.0/eigen_S.min();
   double b_max = -1.0/eigen_Z.min();

   //a_max is de waarde tot waar ik zal zoeken:
   if(b_max < a_max)
      a_max = b_max;

   double a = 0.0;
   double b = a_max;

   double c = (a + b)/2.0;

   //bissectiemethode om stapgrootte te bepalen:
   while(b - a > 1.0e-5){

      c = (a + b)/2.0;

      if( (center_dev + eigen_S.centerpot(c,eigen_Z,c_S,c_Z) - max_dev) < 0.0 )
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * Seperate SUP into two SUP's, a positive and negative semidefinite part.
 * @param p positive (plus) output part
 * @param m negative (minus) output part
 */
void SUP::sep_pm(SUP &p,SUP &m){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sep_pm(p.tpm(i),m.tpm(i));

#ifdef __G_CON

   SZ_ph->sep_pm(p.phm(),m.phm());

#endif

#ifdef __T1_CON

   SZ_dp->sep_pm(p.dpm(),m.dpm());

#endif

#ifdef __T2_CON

   SZ_pph->sep_pm(p.pphm(),m.pphm());

#endif
   
   li->sep_pm(p.gli(),m.gli());

}

/* vim: set ts=3 sw=3 expandtab :*/
