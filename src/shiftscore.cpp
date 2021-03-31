#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec ShiftScoreFast(arma::mat x, arma::mat y, int xn, int yn){
  x /= xn;
  y /= yn;
  
  arma::mat e = arma::cumsum(x) - arma::cumsum(y);
  double ex = e.max();
  double en = e.min();
  ex = ex > 0 ? ex : 0;
  en = en < 0 ? en : 0;
  arma::mat pos = arma::clamp(e, 0, ex);
  arma::mat neg = arma::clamp(e, en, 0);
  
  
  arma::vec ans(4, arma::fill::zeros);
  ans(0) += arma::accu(e);
  ans(1) += arma::accu(pos);
  ans(2) += arma::accu(neg);
  ans(3) += ans(1) - ans(2);
  
  ans /= (e.n_elem - 1);
  return(ans);
}

// [[Rcpp::export]]
arma::vec ShiftScore(arma::mat x, arma::mat y, int calcP, int nresamp){
  // Inputs are sparse vectors of the same length with counts, unnormalized

  int xn = arma::accu(x);
  int yn = arma::accu(y);
  
  arma::vec ans(6, arma::fill::zeros);
  ans.head(4) += ShiftScoreFast(x, y, xn, yn);
  ans(4) += calcP == 0 ? 1.0 : 0.0;
  ans(5) += calcP == 0 ? 1.0 : 0.0;
  
  if (calcP) {
    arma::vec pp(2, arma::fill::zeros);
    
    arma::mat pxy = (x + y) / (xn + yn);
    int k = pxy.n_elem;
    arma::ivec simx(k);
    arma::ivec simy(k);
    arma::vec sim;

    RNGScope scope;
    for (int i = 0; i < nresamp; i++) {
      rmultinom(xn, pxy.begin(), k, simx.begin());
      rmultinom(yn, pxy.begin(), k, simy.begin());
      arma::mat sx = arma::conv_to<arma::mat>::from(simx);
      arma::mat sy = arma::conv_to<arma::mat>::from(simy);
      
      sim = ShiftScoreFast(sx, sy, xn, yn);
      pp(0) += (abs(sim(0)) >= abs(ans(0)));
      pp(1) += (sim(3) >= ans(3));
    }
    ans(4) += pp(0) / nresamp;
    ans(5) += pp(1) / nresamp;
  }


  return(ans);
}


// [[Rcpp::export]]
arma::mat allTheShiftScores(CharacterVector fhash, arma::uvec dists, arma::vec scores,
                            arma::vec sample, int calcP, int nresamp, int ntests){
  arma::mat out(6, ntests);
  out.zeros();

  // loop to find the sequence starts
  int startix = 0;
  int endix = 0;
  int n = fhash.length();
  int ntest = 0;
  for (int it = 1; it < n; it++) {
    if (ntest > ntests - 1) break;
    if (it < n - 1 && fhash(it-1) == fhash(it)) continue;
    endix = it < n-1 ? it - 1 : it;
    if (startix == endix) {
      ntest++;
      continue;
    }
    if (startix >= n - 2 || endix > n - 1) break;
    arma::vec samps = sample.subvec(startix, endix);
    arma::vec vals = scores.subvec(startix, endix);
    arma::uvec samp1 = arma::find(samps);
    arma::uvec samp0 = arma::find(samps < 0.5);
    arma::uvec dists_sub = dists.subvec(startix, endix);
    int nsp = dists_sub.max() + 1;
    arma::uvec locx = dists_sub.elem(samp1);
    arma::uvec locy = dists_sub.elem(samp0);
    
    arma::mat x(nsp, 1, arma::fill::zeros);
    arma::mat y(nsp, 1, arma::fill::zeros);
    x.elem(locx) += vals.elem(samp1);
    y.elem(locy) += vals.elem(samp0);
    
    out.col(ntest) = ShiftScore(x, y, calcP, nresamp);
    ntest++;
    startix = it;
  }  

  return(out);
}
