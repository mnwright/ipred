/*

  $Id: SdiffKM.c,v 1.2 2003/03/27 13:36:02 hothorn Exp $
  
  SdiffKM: integrated squared difference between survival curve 
           and KM estimator

*/
      
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP SdiffKM(SEXP time, SEXP prob, SEXP args)
{
    SEXP rint;
    double d, p, helpone, helptwo, k;
    double myint = 0.0;
    double hazard, window, tw;
    int i, j, n;

    /*
        check arguments
    */

    if (!isVector(time))
	error("Argument time is not a vector");
    n = LENGTH(time);
    if (REAL(time)[0] != 0.0)
        error("time[1] must be zero");

    if (!isVector(prob))
        error("Argument prob is not a vector");
    if (REAL(prob)[0] > 1.0) 
        error("prob[1] must be less or equal 1");
    if (REAL(prob)[0] < 0.0) 
        error("prob[1] must be greater or equal 0");
    if (LENGTH(prob) != n)
        error("prob and time differ in their length");

    if (!isVector(args) || LENGTH(args) != 2)
        error("Argument args must be vector with two elements");

    hazard = REAL(args)[0];
    window = REAL(args)[1];

    /*
        prepare for return values
    */
    
    PROTECT(rint = allocVector(REALSXP, 1));
    UNPROTECT(1);
    REAL(rint)[0] = 0.0;

    /*
        for all discrete times
    */

    for (i = 0; i < n-1; i++) {

      /* get difference between times */
      d = REAL(time)[i+1] - REAL(time)[i];

      /* estimated survival probability at this time */
      p = REAL(prob)[i];

      /* if the difference is small enough */
      if (d < window) {
        helpone = p - exp(-REAL(time)[i] * hazard);
        helptwo = p - exp(-REAL(time)[i+1] * hazard);

        /* mean of over and under sum */
        myint += 0.5 * d * (helpone*helpone + helptwo*helptwo);

      } else {

        /* split up in smaller pieces */
        k = ftrunc(d/window) + 1;
        tw = d/k;
        for (j = 0; j < k; j++) {
          helpone = p - exp(-(REAL(time)[i] + j*tw)*hazard);
          helptwo = p - exp(-(REAL(time)[i] + (j+1)*tw)*hazard);
          
          /* mean of over and under sum for all small windows */
          myint += 0.5*tw*(helpone*helpone + helptwo*helptwo);
        }
      }
    }
    
    /* ok, get outa here */
    
    REAL(rint)[0] = myint;
    return(rint);
}
