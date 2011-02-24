/**
 * @mainpage 
 * This is an implementation of a boundary point method to solve a semidefinite program with linear inequality constraints:
 * we optimizing the second order density matrix using the P Q G T1 and T2 N-representability conditions.
 * At compile time you can decide which condtions will be active compile with make PQ, PQG, PQGT1, PQGT2 or PQGT=(for all conditions).
 * @author Brecht Verstichel, Ward Poelmans
 * @date 23-02-2011
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>

using std::cout;
using std::endl;
using std::ofstream;

#include "include.h"

/**
 * 
 * In the main the actual program is run.\n 
 */

int main(int argc,char **argv)
{
   cout.precision(10);

   // these are the default values
   int M = 8;//dim sp hilbert space
   int N = 4;//nr of particles
   double U = 1;//onsite interaction strength

   struct option long_options[] =
   {
      {"particles",  required_argument, 0, 'n'},
      {"sites",  required_argument, 0, 'm'},
      {"interaction", required_argument, 0, 'U'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;
   while( (j = getopt_long (argc, argv, "hn:m:U:", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -n, --particles=particles    Set the number of particles\n"
               "    -m, --sites=sites            Set the number of sites\n"
               "    -U, --interaction=U          Set the interaction strength\n"
               "    -h, --help                   Display this help\n"
               "\n";
            return 0;
            break;
         case 'n':
            N = atoi(optarg);
            if( N <= 0)
            {
               std::cerr << "Invalid particle number!" << endl;
               return -1;
            }
            break;
         case 'm':
            M = atoi(optarg);
            if( M <= 0)
            {
               std::cerr << "Invalid particle number!" << endl;
               return -2;
            }
            break;
         case 'U':
            U = atof(optarg);
            break;
      }

   cout << "Starting with M=" << M << " N=" << N << " U=" << U << endl;

   LinIneq::init(M,N,1);

   //hamiltoniaan
   TPM ham(M,N);
   ham.hubbard(0,U);

   TPM ham_copy(ham);

   //only traceless hamiltonian needed in program.
   ham.proj_Tr();

   //primal
   SUP X(M,N);

   //dual
   SUP Z(M,N);

   //Lagrange multiplier
   SUP V(M,N);

   //just dubya
   SUP W(M,N);

   SUP u_0(M,N);

   //little help
   TPM hulp(M,N);

   u_0.tpm(0).unit();

   u_0.fill();

   X = 0.0;
   Z = 0.0;

   //what does this do?
   double sigma = 1.6;

   double tolerance = 1.0e-7;

   double D_conv(1.0),P_conv(1.0),convergence(1.0);

   // mazziotti uses 1.6 for this
   double mazzy = 1.0;

   int iter_dual,iter_primal(0);
   int max_iter = 1;

   while(fabs(convergence) > tolerance || P_conv > tolerance || D_conv > tolerance){

      ++iter_primal;

      D_conv = 1.0;

      iter_dual = 0;

      while(D_conv > tolerance  && iter_dual <= max_iter){

         ++iter_dual;

         //solve system
         SUP B(Z);

         B -= u_0;

         B.daxpy(mazzy/sigma,X);

         TPM b(M,N);

         b.collaps(1,B);

         b.daxpy(-mazzy/sigma,ham);

         hulp.S_L(-1,b);

         //hulp is the matrix containing the gamma_i's
         hulp.proj_Tr();

         //construct W
         W.fill(hulp);

         W += u_0;

         W.daxpy(-1.0/sigma,X);

         //update Z and V with eigenvalue decomposition:
         W.sep_pm(Z,V);

         V.dscal(-sigma);

         //check infeasibility of the primal problem:
         TPM v(M,N);

         v.collaps(1,V);

         v -= ham;

         D_conv = sqrt(v.ddot(v));

     }

      //update primal:
      X = V;

      //check dual feasibility (W is a helping variable now)
      W.fill(hulp);

      W += u_0;

      W -= Z;

      P_conv = sqrt(W.ddot(W));
      
      if(D_conv < P_conv)
         sigma *= 1.01;
      else
         sigma /= 1.01;

      convergence = ham.ddot(Z.tpm(0)) + X.ddot(u_0);
         
         cout << P_conv << "\t" << D_conv << "\t" << sigma << "\t" << convergence << endl;

   }

   cout << endl;
   cout << "Energy: " << ham_copy.ddot(Z.tpm(0)) << endl;
   cout << "pd gap: " << Z.ddot(X) << endl;
   cout << "dual conv: " << D_conv << endl;
   cout << "primal conv: " << P_conv << endl;

   LinIneq::clean();

   return 0;

}

/* vim: set ts=3 sw=3 expandtab :*/
