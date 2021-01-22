// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function calculates total per id
// [[Rcpp::export]]
NumericVector GetTotal(IntegerVector id, NumericVector vec, int nobs, int nid) {
  NumericVector tot(nid);

  for (int i = 0; i < nobs; i++) {
    tot[id[i]]=tot[id[i]]+vec[i];
  }
  return tot;
}
