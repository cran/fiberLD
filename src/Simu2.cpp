// Will require the Boost C++ Library for erfc_inv function.

//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>
#include <stdlib.h> // srand
#include <stdio.h> // NULL
//#include <mex.h> 
#include <math.h>
//#include "matrix.h"
#include <time.h> // rand
#include <boost/math/special_functions/erf.hpp>
#include "fiberLD.h"

using namespace Rcpp;


/*
 * Interface to R for Simu2
 */

RcppExport SEXP Simu2Cpp(SEXP x_, SEXP constant_,SEXP mu_,SEXP sigma_,SEXP r_,SEXP r_koll_, SEXP fx_fib_, SEXP antal_,SEXP rnd_seed_)
{
   try{
      double x = as<double>(x_);
      double constant = as<double>(constant_);
      double mu = as<double>(mu_);
      double sigma = as<double>(sigma_);
      double r = as<double>(r_);
      double r_koll = as<double>(r_koll_);
      double fx_fib = as<double>(fx_fib_);
      unsigned int antal = as<unsigned int>(antal_);
      uint32_t rnd_seed = as<uint32_t>(rnd_seed_); 
 
       Rcpp::NumericVector Rout2 = Rcpp::wrap( sim2( x, constant, mu, sigma, r, r_koll, fx_fib, antal, rnd_seed) );

   return Rout2; 
      
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
    return wrap(NA_REAL); 
}



/* 
  *  Internal C++ function for Simu2Cpp
*/


NumericVector sim2(double x, double constant,double mu,double sigma,double r,double r_koll,double fx_fib, unsigned int antal,uint32_t rnd_seed)
{
    double U,hel_sann_fib,w,trial,W,s1,s2;
//    int j,antal,rnd_seed,koll;
    unsigned int i,koll;
    

//do something
    s1=0;
    s2=0;
   // srand(rnd_seed);
    set_seed(rnd_seed);
    
    
    hel_sann_fib=1/(x*sigma*sqrt(2*M_PI))*exp(-(log(x)-mu)*(log(x)-mu)/(2*sigma*sigma))*r_koll/fx_fib;
    
    for (i=0; i<antal; i++) {
        U=Rcpp::runif(1)[0]; // U=((double) rand() / (RAND_MAX));
        if(U<=hel_sann_fib)
        {
            s1=s1+log(x);
            s2=s2+(log(x))*(log(x));
        }
        else{
            koll = 0;
            
            while(koll==0)
            {
                U=Rcpp::runif(1)[0]; //U=((double) rand() / (RAND_MAX));
                W=Rcpp::runif(1)[0];// W=((double) rand() / (RAND_MAX));
                trial=exp(-sqrt(2)*boost::math::erfc_inv(2*(U+(1-U)*0.5*erfc(-(log(x)-mu)/(sqrt(2)*sigma))))*sigma+mu);
                if ((W*constant) <= ((8*r*r-3*x*x+trial*x)/((M_PI*r*r+2*r*trial)*sqrt(4*r*r-x*x))))
                {
                    w=trial;
                    koll=1;
                }
            }
                       
            s1=s1+log(w);
            s2=s2+(log(w))*(log(w));
        }
    }
    NumericVector y2 = NumericVector::create(s1, s2);

    return y2;
}
