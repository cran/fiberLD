/* main method routines */
#include <Rinternals.h>
#include <R.h>
//#include <Rcpp.h>
#include <RcppArmadillo.h>
/* Rconfig.h sometimes doesn't define SUPPORT_OPENMP although
   support is available (e.g. on Windows). Doesn't quite match 
   documentation in `Writing R extensions', but is apparently 
   intentional. However, most compilers with openMP support supply 
   a pre-defined compiler macro _OPENMP. So... */
/* #if (!defined SUPPORT_OPENMP && defined _OPENMP)
  ...update: Rconfig.h:SUPPORT_OPENMP deprecated from R 2.3.2 */
#if defined _OPENMP
#define SUPPORT_OPENMP 1 
#endif

// For safe memory handling from R...
#define CALLOC R_chk_calloc
#define FREE R_chk_free
// Can reset to check for memory errors...
//#define CALLOC calloc
//#define FREE free
Rcpp::NumericVector sim1(double x, double constant,double mu,double sigma,double r,double r_koll,double fx_fin,double zvec, unsigned int antal,uint32_t rnd_seed);

Rcpp::NumericVector sim2(double x, double constant, double mu, double sigma, double r, double r_koll, double fx_fib, unsigned int antal, uint32_t rnd_seed);
	
void set_seed(uint32_t seed);



